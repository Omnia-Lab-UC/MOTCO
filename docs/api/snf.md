# SNF

Similarity Network Fusion and spectral embedding.

## Interpreting trajectory geometry in SNF space

SNF integrates omics layers by constructing sample-affinity graphs, diffusing
information across those graphs, and embedding the fused network. Its natural
geometry is therefore based on neighbourhoods, connectivity, and diffusion.
It does **not** guarantee that feature-space path lengths, vector directions,
or angles are preserved in the spectral coordinates.

MOTCO can use an SNF spectral embedding as its outcome matrix, but the current
trajectory statistics then describe the chosen embedding rather than the
original molecular feature space:

- `delta` is a difference in Euclidean path length in the spectral coordinates;
  it is not automatically a difference in molecular effect magnitude.
- `angle` compares coordinate directions whose relationship to feature-space
  directions is nonlinear and can depend on the fitted graph and eigenspace.
- `shape` compares configurations in the embedding and may respond to changes
  in neighbourhood structure that do not correspond to an isolated molecular
  trajectory bend.

Consequently, Euclidean MOTCO results from SNF should be treated as exploratory
unless simulations demonstrate that the fitted SNF pipeline recovers the
specific geometry of interest. Fit one pooled SNF graph when comparing groups;
separately fitted spectral embeddings do not provide an intrinsically aligned
coordinate system.

For an SNF-native trajectory analysis, prefer quantities defined directly on
the common fused graph, such as:

- diffusion or resistance distance between group-stage cells;
- total diffusion path length and stage-specific transition distances;
- divergence between group-stage transition profiles;
- neighbourhood overlap or graph conductance across successive stages; and
- differences between normalized stage-by-stage diffusion-distance matrices,
  as a scale-free comparison of trajectory configuration.

These are potential methods for a future graph-native test family; MOTCO does
not currently implement them. For feature-space magnitude and orientation,
linear representations such as PCA or PLS are easier to interpret because
their loadings provide an explicit linear map back to the measured features.
A linear space still sees only what lies inside its span: the production
stage-supervised PLS space sees a group difference only through its component
inside the span of the pooled stage means, so a rotation of the stage program
onto other features can vanish from it. See
[Measurement visibility](../measurement-visibility.md).

## Functions

::: motco.stats.snf.get_affinity_matrix

::: motco.stats.snf.SNF

::: motco.stats.snf.get_spectral
