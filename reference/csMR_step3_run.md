# Run csMR Step 3

Run the csMR Snakemake workflow with `R_HOME` and `R_LIBS` explicitly
set from the target csMR conda environment.

## Usage

``` r
csMR_step3_run(
  repo_dir = "~/Project/software/csMR",
  config_file = "./output/csMR/step2_config/config.yml",
  jobs = NULL,
  forcerun = NULL,
  work_flow.snakefile = NULL,
  env_name = "csMR",
  conda = Sys.which("conda"),
  dry_run = FALSE,
  log_file = NULL
)
```

## Arguments

- repo_dir:

  csMR repository directory.

- config_file:

  Config file path.

- jobs:

  Number of Snakemake jobs. If `NULL`, run with bare `-j` and let
  Snakemake use all available cores.

- forcerun:

  Optional Snakemake targets/rules to force-run. Default `NULL` means no
  `--forcerun` is added.

- work_flow.snakefile:

  Snakemake workflow file path. If `NULL`, use
  `file.path(repo_dir, "work_flow.snakefile")`.

- env_name:

  Conda env name.

- conda:

  Path to conda.

- dry_run:

  Whether to run `--dry-run`.

- log_file:

  Optional log file path.

## Value

A list with `command` (copyable shell command), `status`, `log_file`,
`r_lib`, and `r_home`. Output streams to console in real time and is
also saved to `log_file`.

## Details

The first four arguments (`repo_dir`, `config_file`, `jobs`, `forcerun`)
are the most commonly adjusted by users in routine runs.

## Examples

``` r
if (FALSE) { # \dontrun{
csMR_step3_run(repo_dir = "~/Project/software/csMR", jobs = NULL, dry_run = TRUE)
} # }
```
