# Build csMR config.yml

Write an official-style csMR `config.yml` for Step 3.

## Usage

``` r
csMR_step2_config.yml(
  repo_dir = "~/Project/software/csMR",
  config.yml_to = "./output/csMR/step2_config/config.yml",
  BASE_OUTPUT_DIR = "./output/csMR/step3_run",
  qtl_input_dir = "./output/csMR/step1_data_preparation/sc_qtl_dir1",
  exposure_ma = "./output/csMR/step1_data_preparation/gwas_exp_p1.ma",
  exposure_id = "exp_p1",
  exposure_type = "cc",
  exposure_ma_dir = NULL,
  outcome_ma_dir = "./output/csMR/step1_data_preparation/outcome",
  GWAS_REFERENCE_GENOTYPE = NULL,
  eQTL_REFERENCE_GENOTYPE = NULL,
  duplicated_snp_path = "None",
  coloc_window_size_bp = 100000L,
  coloc_coverages = 0.9,
  coloc_threads = 5L,
  coloc_cutoff = 0.8
)
```

## Arguments

- repo_dir:

  csMR repository directory.

- config.yml_to:

  Output config file path.

- BASE_OUTPUT_DIR:

  Output root used by csMR Step 3.

- qtl_input_dir:

  Directory of QTL `.ma` files (or any input supported by
  `.csmr_qtl_paths()`).

- exposure_ma:

  Exposure GWAS `.ma` file path.

- exposure_id:

  Exposure id in `GWAS_SUMSTATS`.

- exposure_type:

  Exposure type, either `"cc"` or `"quant"`.

- exposure_ma_dir:

  Optional directory of multiple exposure `.ma` files. If provided,
  `exposure_ma` is ignored and ids are generated from file names.

- outcome_ma_dir:

  Outcome `.ma` directory (used as `OUTCOME_DIR`). Must contain **only**
  `.ma` files — Snakemake's `os.listdir()` picks up every item (incl.
  `.DS_Store`, subdirectories) as an outcome wildcard.

- GWAS_REFERENCE_GENOTYPE:

  GWAS reference genotype directory. If `NULL`, default to
  `repo_dir/data/reference_genome_1000G_EUR`.

- eQTL_REFERENCE_GENOTYPE:

  eQTL reference genotype directory. If `NULL`, default to
  `repo_dir/data/reference_genome_1000G_EUR`.

- duplicated_snp_path:

  Path to duplicated SNP file for the reference genotype. Use `"None"`
  if not applicable (official csMR default).

- coloc_window_size_bp:

  Coloc window size.

- coloc_coverages:

  Coloc coverage.

- coloc_threads:

  Coloc threads.

- coloc_cutoff:

  Coloc cutoff.

## Value

Absolute path to the written config file.

## Details

Uppercase arguments such as `BASE_OUTPUT_DIR`,
`GWAS_REFERENCE_GENOTYPE`, and `eQTL_REFERENCE_GENOTYPE` intentionally
mirror official csMR `config.yml` keys to keep mapping explicit and
reduce confusion when debugging Step 3.

## Examples

``` r
if (FALSE) { # \dontrun{
out_cfg <- csMR_step2_config.yml(
  repo_dir = "~/Project/software/csMR",
  config.yml_to = "./output/csMR/step2_config/config.yml",
  BASE_OUTPUT_DIR = "./output/csMR/step3_run", # where to store the final csMR output
  qtl_input_dir = "./output/csMR/step1_data_preparation/sc_qtl_dir1",
  exposure_ma = "./output/csMR/step1_data_preparation/gwas_exp_p1.ma",
  exposure_id = "exp_p1",
  exposure_type = "cc",
  outcome_ma_dir = "./output/csMR/step1_data_preparation/outcome"
)
# Example console messages:
# i [23:36:24] Writing csMR config.yml ...
# v [23:36:24] csMR step2 config written: ./output/csMR/step2_config/config.yml
out_cfg
# [1] "/Users/leoarrow/Project/iridocyclitis/output/csMR/step2_config/config.yml"
} # }
```
