# Adjust SMR Results with FDR and Bonferroni Corrections

It applies FDR/Bonferroni corrections for a single SMR result. The
corrected results are saved in an `fdr` subdirectory within the output
directory.

## Usage

``` r
leo_smr_adjust(
  smr_result_path,
  writePath = "",
  out_dir = "",
  QTL_type = "",
  Source = "",
  Tissue = "",
  Outcome_name = "",
  add_info_cols = T,
  drop_non_heidi = T,
  hla = F,
  write_out = T
)
```

## Arguments

- smr_result_path:

  Character. The path containing SMR result file. Files should have the
  `.smr` extension.

- writePath:

  Character. The path to write the adjusted results.

- out_dir:

  Character. The output directory where adjusted files will be saved. If
  not specified, defaults to create a dir named `fdr` within dir-name
  (smr_result_path).

- QTL_type:

  Character. The type of the QTL for SMR analysis.

- Source:

  Character. The name of the Source.

- Tissue:

  Character. The name of the Tissue.

- Outcome_name:

  Character. The name of the Outcome.

- add_info_cols:

  Logical. Whether to add information columns for: QTL_type, Source,
  Tissue and Outcome. Default is `TRUE`.

- drop_non_heidi:

  Logical. Whether to drop Probes without HEIDI test information.
  Default is `TRUE`.

- hla:

  Logical. Whether to pre-exclude the HLA region probes. Default is
  `False`.

- write_out:

  Logical. Whether to write the adjusted results to a file. Default is
  `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
leo_smr_adjust(
  "~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined/chr_combine_sQTL_Adipose_Subcutaneous@iri3.smr",
  writePath = "", out_dir = "~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49"
)
leo_smr_adjust(
  paste0("~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined/",
         "chr_combine_sQTL_Adipose_Subcutaneous@iri3.smr"),
  writePath = "./haha.fdr", out_dir = ""
)
} # }
```
