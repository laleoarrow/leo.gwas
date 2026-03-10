# Configure Conda Environment for csMR

`csMR_env()` prepares the official csMR environment. It clones the
official repository if missing, tries `envs/envpy3.yml` first, falls
back to a temporary compatibility env file only when the official file
cannot solve, installs required R packages into the csMR conda library,
checks PLINK 1.90, and returns the runtime paths required by the csMR
README.

## Usage

``` r
csMR_env(
  repo_dir = "~/Project/software/csMR",
  repo_url = "https://github.com/rhhao/csMR.git",
  env_name = "csMR",
  conda = Sys.which("conda"),
  mamba = Sys.which("mamba"),
  plink_path = Sys.which("plink"),
  ref = "main",
  overwrite = FALSE
)
```

## Arguments

- repo_dir:

  Path to the csMR repository.

- repo_url:

  Official csMR GitHub URL.

- env_name:

  Conda environment name.

- conda:

  Path to conda.

- mamba:

  Path to mamba.

- plink_path:

  Optional PLINK 1.90 path. If empty, `PATH` and
  `plinkbinr::get_plink_exe()` are tried.

- ref:

  Git branch or tag.

- overwrite:

  Whether to rebuild an existing env.

## Value

A list with `repo_dir`, `env_name`, `env_prefix`, `env_file`, `r_home`,
`r_library`, and `plink`.

## Examples

``` r
if (FALSE) { # \dontrun{
cfg <- csMR_env()
cfg$r_library
} # }
```
