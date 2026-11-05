"""Statistical utilities for MOTCO.

Modules
-------
pls
    Partial Least Squares (PLS-DA) utilities.
snf
    Similarity Network Fusion and spectral embedding.
design
    Design matrix and LS-means construction.
trajectory
    Trajectory estimation and geometric metrics.
permutation
    RRPP permutation infrastructure.
interaction
    Group × stage interaction test split into a linear span and its orthogonal complement.
attribution
    Feature attribution for two-group trajectories in a shared PLS space.
"""

from .attribution import (  # noqa: F401
    AggregateEffect,
    AttributionConfig,
    AttributionError,
    AttributionResult,
    BootstrapSummary,
    ComponentAttribution,
    FeatureEffect,
    InterpretationMetadata,
    OrientationAttributionResult,
    PrincipalOrientationAttribution,
    TransitionAttribution,
    analyze_orientation_attribution,
    attribution_frames,
    run_orientation_attribution,
    write_attribution_outputs,
)
from .design import build_ls_means, center_matrix, get_model_matrix  # noqa: F401
from .interaction import InteractionSplitPart, InteractionSplitResult, interaction_split_test  # noqa: F401
from .permutation import RRPP  # noqa: F401
from .pls import calculate_vips, fit_plsda_model, fit_plsda_transform, plsda_doubleCV  # noqa: F401
from .snf import SNF, get_affinity_matrix, get_spectral  # noqa: F401
from .trajectory import (  # noqa: F401
    configuration_spectra,
    configuration_spectrum,
    estimate_betas,
    estimate_difference,
    get_observed_vectors,
    pair_difference,
    principal_orientation,
)

__all__ = [
    # pls
    "plsda_doubleCV",
    "calculate_vips",
    "fit_plsda_model",
    "fit_plsda_transform",
    # design
    "build_ls_means",
    "center_matrix",
    "get_model_matrix",
    # trajectory
    "configuration_spectra",
    "configuration_spectrum",
    "estimate_betas",
    "estimate_difference",
    "get_observed_vectors",
    "pair_difference",
    "principal_orientation",
    # attribution
    "AggregateEffect",
    "AttributionConfig",
    "AttributionError",
    "AttributionResult",
    "BootstrapSummary",
    "ComponentAttribution",
    "FeatureEffect",
    "InterpretationMetadata",
    "OrientationAttributionResult",
    "PrincipalOrientationAttribution",
    "TransitionAttribution",
    "analyze_orientation_attribution",
    "attribution_frames",
    "run_orientation_attribution",
    "write_attribution_outputs",
    # permutation
    "RRPP",
    # interaction
    "InteractionSplitPart",
    "InteractionSplitResult",
    "interaction_split_test",
    # snf
    "SNF",
    "get_affinity_matrix",
    "get_spectral",
]
