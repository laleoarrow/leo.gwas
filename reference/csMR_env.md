# Configure Conda Environment for csMR

`csMR_env()` automates environment setup for the official csMR pipeline:
it clones/updates the csMR repository, creates (or updates) a conda
environment from `envs/envpy3.yml`, installs missing R dependencies used
by csMR scripts, and runs basic command/package checks.

## Usage

``` r
csMR_env(
  repo_dir = file.path(path.expand("~"), "Project", "software"),
  env_name = "csMR",
  conda = Sys.which("conda"),
  repo_url = "https://github.com/rhhao/csMR.git",
  ref = "main",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- repo_dir:

  Directory where the csMR repository is stored. You can pass either the
  csMR repository path itself or its parent directory (for example
  `~/Project/software`, resolved to `~/Project/software/csMR`).

- env_name:

  Name of the conda environment.

- conda:

  Path to the conda executable. Defaults to `Sys.which("conda")`.

- repo_url:

  GitHub URL of the csMR repository.

- ref:

  Git branch/tag to checkout for csMR source code.

- overwrite:

  Whether to remove an existing environment before recreating.

- verbose:

  Whether to print setup logs.

## Value

A list containing `repo_dir`, `env_name`, `env_prefix`, `env_file`,
`r_library`, and two convenience shell snippets (`export_r_libs`,
`run_snakemake_example`).

## Details

The method follows the csMR paper workflow (Nat Commun 2024, 15:4890):
(1) locus-level precomputation from GWAS and single-cell eQTL summary
data, (2) SuSiE-based fine-mapping and colocalization, and (3)
cell-type-specific MR with sensitivity analyses. This function only
prepares software/runtime dependencies; it does not execute the full
csMR analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# First-time setup
cfg <- csMR_env(repo_dir = "~/Project/software", env_name = "csMR")
cfg$run_snakemake_example

# Fast health-check for an existing environment
csMR_env(repo_dir = "~/Project/software", env_name = "csMR")
} # }
```
