# Prepare csMR Step1 Input

Convert GWAS or QTL summary statistics to the csMR-required `.ma` format
using explicit column mappings. This function only performs thin
formatting and basic QC; it does not guess genome build or map non-rsID
SNPs.

## Usage

``` r
csMR_step1_prep(
  data,
  type = c("gwas", "qtl"),
  output,
  SNP_col = "SNP",
  GENE_col = "GENE",
  A1_col = "A1",
  A2_col = "A2",
  MAF_col = "MAF",
  BETA_col = "BETA",
  SE_col = "SE",
  P_col = "P",
  N_col = "N",
  n = NULL
)
```

## Arguments

- data:

  GWAS/QTL data.frame, file path, file vector, or a QTL directory.

- type:

  `"gwas"` or `"qtl"`.

- output:

  Output `.ma` file path for GWAS or single-table QTL, or output
  directory for multi-file QTL input.

- SNP_col:

  Input SNP column name.

- GENE_col:

  Input gene column name for QTL.

- A1_col:

  Input effect allele column name.

- A2_col:

  Input other allele column name.

- MAF_col:

  Input MAF column name.

- BETA_col:

  Input beta column name.

- SE_col:

  Input SE column name.

- P_col:

  Input P column name.

- N_col:

  Input sample size column name.

- n:

  Optional fixed sample size used to fill missing `N` in GWAS.

## Value

For GWAS, a list with `output`, `n_input`, `n_output`, and
`n_missing_filled`. For QTL, a manifest data.frame with `id`, `input`,
`output`, `n_input`, and `n_output`.

## Examples

``` r
if (FALSE) { # \dontrun{
csMR_step1_prep(
  data = "~/Project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv",
  type = "gwas",
  output = "~/Project/iridocyclitis/output/csMR/step1/exposure.ma",
  SNP_col = "rsID",
  A1_col = "effect_allele",
  A2_col = "other_allele",
  MAF_col = "EAF",
  BETA_col = "beta",
  SE_col = "se",
  P_col = "pval",
  N_col = "N"
)
} # }
```
