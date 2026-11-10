# Reference fixture for motco.stats.sparse_pls.SparsePLSDA (regularized-stage-pls).
#
# Fits mixOmics::splsda (regression mode, no per-feature rescaling) on a seeded
# three-class Gaussian X and writes the inputs and the fitted weights and
# scores as CSV. Run once with an R that has mixOmics, e.g. the flake's dev
# shell (`nix develop`), from the repository root:
#
#     Rscript tests/data/reference/splsda_fixture.R build/splsda_fixture
#     uv run python -m tests.splsda_fixture_capture build/splsda_fixture
#
# The capture script packs tests/data/splsda_reference_fixture.npz with this
# script's sha256 and the R and mixOmics versions; tests never need R.

suppressPackageStartupMessages(library(mixOmics))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript splsda_fixture.R <out_dir>")
out_dir <- args[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20261110)
n_per_class <- 15
n_features <- 60
keep <- 8
ncomp <- 3

classes <- rep(c("a", "b", "c"), each = n_per_class)
X <- matrix(rnorm(length(classes) * n_features), nrow = length(classes))
# Class signal on the first twelve features, in two directions.
X[, 1:6] <- X[, 1:6] + 1.2 * (classes == "b")
X[, 7:12] <- X[, 7:12] + 0.9 * (classes == "c")
colnames(X) <- sprintf("f%02d", seq_len(n_features))

fit <- splsda(
  X, factor(classes), ncomp = ncomp, keepX = rep(keep, ncomp),
  scale = FALSE, tol = 1e-20, max.iter = 10000
)

write.csv(X, file.path(out_dir, "X.csv"), row.names = FALSE)
write.csv(data.frame(class = classes), file.path(out_dir, "classes.csv"), row.names = FALSE)
write.csv(fit$ind.mat, file.path(out_dir, "Y.csv"), row.names = FALSE)
write.csv(fit$loadings$X, file.path(out_dir, "weights.csv"), row.names = FALSE)
write.csv(fit$variates$X, file.path(out_dir, "scores.csv"), row.names = FALSE)
writeLines(
  c(
    paste0("r_version=", R.version.string),
    paste0("mixomics_version=", as.character(packageVersion("mixOmics"))),
    paste0("keep=", keep),
    paste0("ncomp=", ncomp),
    paste0("iterations=", paste(fit$iter, collapse = ","))
  ),
  file.path(out_dir, "provenance.txt")
)
