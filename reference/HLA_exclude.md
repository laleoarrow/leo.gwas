# Exclude HLA region from genomic summary data

This function filters out entries within the HLA region on a specified
chromosome and position range. The default HLA region is set to
chromosome 6, between 25Mb and 34Mb. Custom chromosome and position
bounds can be specified.

## Usage

``` r
HLA_exclude(
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

A data frame excluding rows that fall within the specified HLA region.

## Examples

``` r
example_data <- data.frame(
  SNP = c("rs1", "rs2", "rs3", "rs4"),
  CHR = c(6, 6, 6, 7),
  POS = c(26000000, 33000000, 35000000, 29000000)
)
result_data <- HLA_exclude(example_data, chromosome_col="CHR", position_col="POS")
print(result_data)
#>   SNP CHR     POS
#> 1 rs3   6 3.5e+07
#> 2 rs4   7 2.9e+07
```
