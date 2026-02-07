# Extract HLA region from genomic summary data

This function extracts entries within the HLA region on a specified
chromosome and position range. The default HLA region is chromosome 6,
between 25Mb and 34Mb. Users may provide custom chromosome and region
boundaries.

## Usage

``` r
HLA_get(
  data,
  chromosome_col = "CHR",
  position_col = "BP",
  lower_bound = 2.5e+07,
  upper_bound = 3.4e+07
)
```

## Arguments

- data:

  A data frame containing genomic data.

- chromosome_col:

  The name of the column representing chromosome numbers (default is
  CHR).

- position_col:

  The name of the column representing genomic positions (default is BP).

- lower_bound:

  The lower boundary of the HLA region in base pairs (default is 25e6).

- upper_bound:

  The upper boundary of the HLA region in base pairs (default is 34e6).

## Value

A data frame including only rows within the specified HLA region.

## Examples

``` r
example_data <- data.frame(
  SNP = c("rs1", "rs2", "rs3", "rs4"),
  CHR = c(6, 6, 6, 7),
  POS = c(26000000, 33000000, 35000000, 29000000)
)
hla_data <- HLA_get(example_data, chromosome_col="CHR", position_col="POS")
print(hla_data)
#>   SNP CHR     POS
#> 1 rs1   6 2.6e+07
#> 2 rs2   6 3.3e+07
```
