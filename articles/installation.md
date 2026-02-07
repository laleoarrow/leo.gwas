# Installation Tutorial

## Aim

Install `leo.gwas` with the most robust path first, then use fallback
commands when specific dependencies fail.

## Recommended path (2026-02)

``` r
options(timeout = 1000000)

# 1) Install leo.basic first
pak::pkg_install("github::laleoarrow/leo.basic")

# 2) Install heavy/optional dependencies through leo.basic helper
# If your helper function name is install.deps, use this:
leo.basic::install.deps("leo.gwas")

# If your local leo.basic uses install_deps naming, use this instead:
# leo.basic::install_deps("leo.gwas")

# 3) Install leo.gwas
pak::pkg_install("github::laleoarrow/leo.gwas")
# or
# devtools::install_github("laleoarrow/leo.gwas")
```

## Fallback commands for common failures

If dependency resolution is still unlucky, run the following manually
and retry `leo.gwas` installation.

``` r
# MR-related dependencies
remotes::install_github("MRCIEU/TwoSampleMR")

# PLINK helper (optional; needed for local clumping workflows)
# This package may occasionally fail on CI/source unpack in some environments.
devtools::install_github("explodecomputer/plinkbinr")

# Bioconductor dependency example
BiocManager::install("minfi")

# CatBoost (official binary-installation page)
# https://catboost.ai/docs/en/installation/r-installation-binary-installation

# Recommended binary install command as of 2026-02 (macOS universal2)
pak::pak("url::https://github.com/catboost/catboost/releases/download/v1.2.8/catboost-R-darwin-universal2-1.2.8.tgz")

# Equivalent remotes command example
remotes::install_url(
  "https://github.com/catboost/catboost/releases/download/v1.2.8/catboost-R-darwin-universal2-1.2.8.tgz",
  INSTALL_opts = c("--no-multiarch", "--no-test-load", "--no-staged-install")
)

# ggrastr fallback
devtools::install_github("VPetukhov/ggrastr", build_vignettes = FALSE)
```

## Verify installation

``` r
library(leo.gwas)
packageVersion("leo.gwas")
```

If [`library(leo.gwas)`](https://laleoarrow.github.io/leo.gwas/) and
`packageVersion("leo.gwas")` run successfully, installation is complete.

## Notes

- `catboost`, `caret`, `plinkbinr`, and some genomics packages are
  optional for specific functions.
- If a function reports missing package errors, install only that
  package and rerun.
