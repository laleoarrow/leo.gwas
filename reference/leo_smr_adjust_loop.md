# Batch Adjust SMR Results with FDR and Bonferroni Corrections

This function applies FDR/Bonferroni corrections to all SMR results
within a specified directory.

## Usage

``` r
leo_smr_adjust_loop(
  dir,
  out_dir = "",
  pattern = "\\.smr$",
  QTL_type,
  Source,
  ...
)
```

## Arguments

- dir:

  Character. The directory containing SMR result files. Files should
  have the `.smr` extension and follow the naming convention
  `exposure@outcome.smr`.

- out_dir:

  Character. The output directory where the adjusted SMR files will be
  saved. If not specified, defaults to creating a `fdr` directory within
  the input directory.

- pattern:

  Character. A regular expression pattern to match SMR result files.
  Default is `"\.smr$"` (i.e., files ending with `.smr`).

- QTL_type:

  Character. Type of QTL (e.g., "eQTL", "sQTL").

- Source:

  Character. Source of the data (e.g., "GTEx", "Westra").

- ...:

  Additional arguments to be passed to `leo_smr_adjust`.

## Value

NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
leo_smr_adjust_loop(dir     = "./output/smr/sqtl/GTEx49/chr_combined",
                    out_dir = "./output/smr/sqtl/GTEx49/fdr")
} # }
```
