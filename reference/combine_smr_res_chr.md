# Combine SMR Results for All Chromosomes

This function combines SMR (Summary-data-based Mendelian Randomization)
result files across all chromosomes for each unique exposure and outcome
pair (Yes, we can deal with multiple exposure in one dir for say, the 49
sqtl from GTEx). The combined results are saved to a specified output
directory.

## Usage

``` r
combine_smr_res_chr(dir, out_dir = "")
```

## Arguments

- dir:

  Character. The directory containing SMR result files. Files should
  follow the naming convention `exposure_chrX@outcome.smr`, where `X`
  represents the chromosome number.

- out_dir:

  Character. The output directory where combined SMR files will be
  saved. If not specified, defaults to a subdirectory named
  `chr_combined` within `dir`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Combine SMR results in the "data/smr_results" directory and save to default output directory
combine_smr_res_chr(dir = "data/smr_results")

# Combine SMR results and specify a custom output directory
combine_smr_res_chr(dir = "data/smr_results", out_dir = "data/combined_results")
} # }
```
