"""numpy-native multi-omic generator reproducing InterSIM's mean-shift model.

This replaces the per-replicate R subprocess. InterSIM's generative model
(read from ``body(InterSIM)``) is, per omic and per cluster ``i``::

    effect_i = base + delta * v_i          # v_i in {0,1}^p is the differential indicator
    X_i      ~ MVN(effect_i, Sigma)        # methylation then passes through rev.logit

with cross-omic coupling: differential genes are derived from differential CpGs
and differential proteins from differential genes, and expression/protein means
blend a per-feature cross-omic correlation term. All reference objects come from
the cached :class:`~motco.simulations.reference.IntersimReference`; no R is used
at runtime.

The trajectory layer treats each *cell* (a group x stage combination) as one
"cluster" column, so the same core handles InterSIM-style clustering and the
feature-surgery trajectory modes uniformly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from motco.simulations.reference import IntersimReference, load_reference

#: Omic layer names in canonical order.
OMIC_LAYERS = ("methylation", "expression", "proteomics")


class GeneratorError(ValueError):
    """Raised for invalid generator inputs."""


@dataclass(frozen=True)
class GeneratedOmics:
    """Aligned omics matrices plus the indicator truth used to generate them.

    ``methylation`` / ``expression`` / ``proteomics`` are ``(n_samples, n_feat)``
    matrices sharing the row order given by ``cell_ids`` (the index of the cell
    each sample was drawn from). ``indicators_*`` are ``(n_feat, n_cell)`` {0,1}
    matrices recording the differential-feature indicator for each cell.
    """

    methylation: np.ndarray
    expression: np.ndarray
    proteomics: np.ndarray
    cell_ids: np.ndarray
    indicators_methyl: np.ndarray
    indicators_expr: np.ndarray
    indicators_protein: np.ndarray


def omic_population_means(
    *,
    indicators_methyl: np.ndarray,
    indicators_expr: np.ndarray,
    indicators_protein: np.ndarray,
    delta_methyl: float = 2.0,
    delta_expr: float = 2.0,
    delta_protein: float = 2.0,
    reference: IntersimReference | None = None,
    base_shift_methyl: np.ndarray | None = None,
    base_shift_expr: np.ndarray | None = None,
    base_shift_protein: np.ndarray | None = None,
    cell_shift_methyl: np.ndarray | None = None,
    cell_shift_expr: np.ndarray | None = None,
    cell_shift_protein: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Return exact cell means in the units supplied to integration.

    Arrays have shape ``(n_cells, n_features)``. Methylation means are in
    pre-``rev_logit`` M-value units; expression and proteomics retain their
    native generated units. ``cell_shift_*`` are optional per-cell mean shifts
    shaped like the indicators, ``(n_features, n_cells)``, in the same units
    (methylation in M-values); an absent shift leaves the means untouched.
    """

    ref = reference if reference is not None else load_reference()
    rho_me = np.zeros_like(ref.rho_methyl_expr) if delta_expr == 0 else ref.rho_methyl_expr
    base_expr = rho_me * ref.methyl_gene_level_mean + np.sqrt(1 - rho_me**2) * ref.mean_expr
    rho_ep = np.zeros_like(ref.rho_expr_protein) if delta_protein == 0 else ref.rho_expr_protein
    base_protein = rho_ep * ref.mean_expr_with_mapped_protein + np.sqrt(1 - rho_ep**2) * ref.mean_protein
    return {
        "methylation": _cell_means(
            ref.mean_M, delta_methyl, indicators_methyl, base_shift_methyl, cell_shift_methyl
        ),
        "expression": _cell_means(base_expr, delta_expr, indicators_expr, base_shift_expr, cell_shift_expr),
        "proteomics": _cell_means(
            base_protein, delta_protein, indicators_protein, base_shift_protein, cell_shift_protein
        ),
    }


def rev_logit(x: np.ndarray) -> np.ndarray:
    """InterSIM's inverse-logit: ``1 / (1 + exp(-x))``."""

    return 1.0 / (1.0 + np.exp(-x))


def logit(x: np.ndarray, clip: float = 1e-6) -> np.ndarray:
    """Clipped logit (M-value): exact inverse of ``rev_logit``.

    Clips ``x`` to ``[clip, 1 - clip]`` before applying ``ln(x / (1 - x))``
    to guard against divergence when B values are numerically at 0 or 1.
    """

    b = np.clip(x, clip, 1.0 - clip)
    return np.log(b / (1.0 - b))


def bernoulli_indicators(
    rng: np.random.Generator, n_feat: int, n_cell: int, p: float
) -> np.ndarray:
    """Independent Bernoulli(``p``) differential indicator per feature and cell."""

    return (rng.random((n_feat, n_cell)) < p).astype(float)


def markov_indicators(
    rng: np.random.Generator, n_feat: int, n_cell: int, p: float, rho: float = 0.0
) -> np.ndarray:
    """Differential indicators that persist across *ordered* cells with continuity ``rho``.

    Each feature's indicator sequence follows a stationary first-order Markov
    chain along the cell axis (which the trajectory layer orders as stages)::

        x_1 ~ Bern(p)
        P(x_t = 1 | x_{t-1} = 1) = p + rho * (1 - p)      # stay differential
        P(x_t = 1 | x_{t-1} = 0) = p * (1 - rho)          # become differential

    **Stationary-marginal contract:** the chain's stationary distribution is
    Bernoulli(``p``) and the first cell is drawn from it, so *every* cell's
    marginal is Bernoulli(``p``) at every ``rho``. Per-cell differential counts,
    the per-omic effect sizes, and the CpG→gene→protein coupling derivation are
    therefore invariant along the continuity axis; only the *cross-cell*
    correlation changes, as ``corr(x_t, x_s) = rho ** abs(t - s)``.

    ``rho = 0`` collapses both transition probabilities to ``p``, and the single
    uniform block drawn here has the same shape and stream position as
    :func:`bernoulli_indicators`, so the output is *identical* to that function's
    on the same fresh generator and nothing drawn afterwards shifts.
    """

    if not (0.0 <= rho < 1.0):
        raise GeneratorError(f"rho must be in [0, 1); got {rho}.")
    uniforms = rng.random((n_feat, n_cell))
    if n_cell == 0:
        return uniforms.astype(float)
    stay = p + rho * (1.0 - p)
    enter = p * (1.0 - rho)
    indicators = np.empty((n_feat, n_cell), dtype=float)
    previous = uniforms[:, 0] < p
    indicators[:, 0] = previous
    for cell in range(1, n_cell):
        previous = uniforms[:, cell] < np.where(previous, stay, enter)
        indicators[:, cell] = previous
    return indicators


def derive_coupled_indicators(
    indicators_methyl: np.ndarray, reference: IntersimReference
) -> tuple[np.ndarray, np.ndarray]:
    """Derive expression and protein indicators from methylation indicators.

    Mirrors InterSIM's default coupling: a gene is differential when any of its
    mapped CpGs is, and a protein is differential when its mapped gene is. Uses
    the cached incidence matrices, so the result matches InterSIM's
    ``CpG.gene.map.for.DEG`` / ``protein.gene.map.for.DEP`` derivation.
    """

    expr = (reference.incidence_cpg_gene.T @ indicators_methyl > 0).astype(float)
    protein = (reference.incidence_gene_protein.T @ expr > 0).astype(float)
    return expr, protein


def _sample_cells(
    rng: np.random.Generator,
    base: np.ndarray,
    delta: float,
    indicators: np.ndarray,
    cov: np.ndarray,
    cell_sizes: list[int],
    base_shift: np.ndarray | None = None,
) -> np.ndarray:
    """Stack MVN draws for each cell with mean ``base + base_shift + delta * indicator``.

    ``base_shift`` is a constant per-feature offset added to every cell's mean
    (used for the location-only ``translation`` mode, applied in the pre-logit
    M-value space for methylation).
    """

    means = _cell_means(base, delta, indicators, base_shift)
    blocks = []
    for i, n in enumerate(cell_sizes):
        if n == 0:
            continue
        blocks.append(rng.multivariate_normal(means[i], cov, size=n, method="svd"))
    return np.vstack(blocks)


def _sample_population_means(
    rng: np.random.Generator,
    means: np.ndarray,
    cov: np.ndarray,
    cell_sizes: list[int],
) -> np.ndarray:
    """Draw each cell from its corresponding exact population mean."""

    blocks = [
        rng.multivariate_normal(means[i], cov, size=n, method="svd")
        for i, n in enumerate(cell_sizes)
        if n > 0
    ]
    return np.vstack(blocks)


def _cell_means(
    base: np.ndarray,
    delta: float,
    indicators: np.ndarray,
    base_shift: np.ndarray | None = None,
    cell_shift: np.ndarray | None = None,
) -> np.ndarray:
    """Construct exact cell means with cells on rows and features on columns.

    ``cell_shift`` (``(n_features, n_cells)``, like ``indicators``) is added to
    each cell's mean last; ``None`` takes the path used before it existed.
    """

    shifted_base = base if base_shift is None else base + base_shift
    means = shifted_base[None, :] + delta * indicators.T
    if cell_shift is None:
        return means
    shift = np.asarray(cell_shift, dtype=float)
    if shift.shape != indicators.shape:
        raise GeneratorError(f"cell shift must have shape {indicators.shape}; got {shift.shape}")
    return means + shift.T


def _resolve_cov(
    sigma: str | np.ndarray | None, default: np.ndarray
) -> np.ndarray:
    if sigma is None:
        return default
    if isinstance(sigma, str):
        if sigma == "indep":
            return np.diag(np.diag(default))
        raise GeneratorError(f"Unknown covariance option: {sigma!r}")
    return np.asarray(sigma, dtype=float)


def generate_omics(
    *,
    cell_sizes: list[int],
    indicators_methyl: np.ndarray,
    indicators_expr: np.ndarray,
    indicators_protein: np.ndarray,
    delta_methyl: float = 2.0,
    delta_expr: float = 2.0,
    delta_protein: float = 2.0,
    rng: np.random.Generator,
    reference: IntersimReference | None = None,
    sigma_methyl: str | np.ndarray | None = None,
    sigma_expr: str | np.ndarray | None = None,
    sigma_protein: str | np.ndarray | None = None,
    base_shift_methyl: np.ndarray | None = None,
    base_shift_expr: np.ndarray | None = None,
    base_shift_protein: np.ndarray | None = None,
    cell_shift_methyl: np.ndarray | None = None,
    cell_shift_expr: np.ndarray | None = None,
    cell_shift_protein: np.ndarray | None = None,
) -> GeneratedOmics:
    """Generate aligned omics matrices from explicit per-cell indicators.

    Each "cell" (column of the indicator matrices) is sampled as a cluster with
    mean ``base + delta * indicator`` and the reference covariance, exactly as
    InterSIM does. Methylation is shifted in M-value (logit) space and then
    passed through :func:`rev_logit`; expression and protein means carry the
    InterSIM cross-omic blend. ``cell_shift_*`` add an optional per-cell mean
    shift (``(n_features, n_cells)``; methylation in M-values, before
    :func:`rev_logit`) and draw no random numbers, so at a fixed ``rng`` the
    samples move by exactly the shift.
    """

    ref = reference if reference is not None else load_reference()
    n_cell = len(cell_sizes)
    for name, ind, n_feat in (
        ("methylation", indicators_methyl, ref.n_cpg),
        ("expression", indicators_expr, ref.n_gene),
        ("proteomics", indicators_protein, ref.n_protein),
    ):
        if ind.shape != (n_feat, n_cell):
            raise GeneratorError(
                f"{name} indicators must have shape ({n_feat}, {n_cell}); got {ind.shape}"
            )
    if any(n < 0 for n in cell_sizes):
        raise GeneratorError("cell_sizes must be non-negative")

    cov_M = _resolve_cov(sigma_methyl, ref.cov_M)
    cov_expr = _resolve_cov(sigma_expr, ref.cov_expr)
    cov_protein = _resolve_cov(sigma_protein, ref.cov_protein)

    means = omic_population_means(
        indicators_methyl=indicators_methyl,
        indicators_expr=indicators_expr,
        indicators_protein=indicators_protein,
        delta_methyl=delta_methyl,
        delta_expr=delta_expr,
        delta_protein=delta_protein,
        reference=ref,
        base_shift_methyl=base_shift_methyl,
        base_shift_expr=base_shift_expr,
        base_shift_protein=base_shift_protein,
        cell_shift_methyl=cell_shift_methyl,
        cell_shift_expr=cell_shift_expr,
        cell_shift_protein=cell_shift_protein,
    )

    # Methylation: additive shift in M-value space, then inverse-logit.
    methyl_M = _sample_population_means(rng, means["methylation"], cov_M, cell_sizes)
    methylation = rev_logit(methyl_M)

    expression = _sample_population_means(rng, means["expression"], cov_expr, cell_sizes)
    proteomics = _sample_population_means(rng, means["proteomics"], cov_protein, cell_sizes)

    cell_ids = np.concatenate([np.full(n, i) for i, n in enumerate(cell_sizes) if n > 0])

    return GeneratedOmics(
        methylation=methylation,
        expression=expression,
        proteomics=proteomics,
        cell_ids=cell_ids,
        indicators_methyl=indicators_methyl,
        indicators_expr=indicators_expr,
        indicators_protein=indicators_protein,
    )
