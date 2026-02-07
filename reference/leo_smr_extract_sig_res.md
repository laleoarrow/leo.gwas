# Extract significant results from `.all` files

This function scans the `combine_1outcome` folder for `.all` files,
filters rows where `Pass_FDR == "Pass"` or `Pass_Bonferroni == "Pass"`
(depending on `pass_type` argument), and writes the significant subset
to a new file (e.g., `smr_res_YYYYMMDD_iri3.sig.all`).

## Usage

``` r
leo_smr_extract_sig_res(dir, pass_type = c("FDR", "Bonferroni"), out_dir = "")
```

## Arguments

- dir:

  Character. The main `combine_1outcome` folder from
  [combine_smr_res_1outcome](https://laleoarrow.github.io/leo.gwas/reference/combine_smr_res_1outcome.md).

- pass_type:

  Character. Which significance criterion to use: one of
  `c("FDR", "Bonferroni", "both")`. - "FDR": keep rows where
  `Pass_FDR == "Pass"` - "Bonferroni": keep rows where
  `Pass_Bonferroni == "Pass"` - "both": keep rows where either FDR or
  Bonferroni is "Pass"

- out_dir:

  Character. Where to write significant results. Default
  `dir/combine_1outcome/sig`.

## Value

NULL. Writes `.sig.all` files to disk.

## Examples

``` r
if (FALSE) { # \dontrun{
# For the directory "/Users/leoarrow/project/iridocyclitis/output/smr-t2d",
# after running combine_smr_res_1outcome, we have .all files in
# "smr-t2d/combine_1outcome".

leo_smr_extract_sig_res(
  dir       = "/Users/leoarrow/project/iridocyclitis/output/smr-t2d",
  pass_type = "FDR"
)
} # }
```
