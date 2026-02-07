# Filter Chromosomes Based on SNP P-value Threshold

This function filters chromosomes out if no SNP within the chromosome
meets threshold. That is, if a chromosome has at least one SNP that
meets the threshold, all SNPs on that chromosome are retained.

## Usage

``` r
filter_chr_basedonSNP_p(
  df,
  chr_col = "CHR",
  snp_col = "Variant_ID",
  p_val_col = "nominal_P_value",
  threshold = 0.00157
)
```

## Arguments

- df:

  A data frame containing SNP data.

- chr_col:

  Character string specifying the name of the chromosome column. Default
  is `"CHR"`.

- snp_col:

  Character string specifying the name of the SNP identifier column.
  Default is `"Variant_ID"`.

- p_val_col:

  Character string specifying the name of the p-value column. Default is
  `"nominal_P_value"`.

- threshold:

  Numeric value specifying the p-value threshold. Default is `1.57e-3`
  (Thresh hold for the HEIDI test).

## Value

A filtered data frame .

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
eqtl_data <- data.frame(
CHR = c("1", "1", "2", "2", "3"),
Variant_ID = c("rs1", "rs2", "rs3", "rs4", "rs5"),
nominal_P_value = c(1e-9, 0.05, 0.2, 1e-7, 0.3)
);eqtl_data
#>   CHR Variant_ID nominal_P_value
#> 1   1        rs1           1e-09
#> 2   1        rs2           5e-02
#> 3   2        rs3           2e-01
#> 4   2        rs4           1e-07
#> 5   3        rs5           3e-01
df_filtered <- filter_chr_basedonSNP_p(
  df = eqtl_data,
  chr_col = "CHR",
  snp_col = "Variant_ID",
  p_val_col = "nominal_P_value",
  threshold = 5e-8
);df_filtered
#> # A tibble: 2 × 3
#>   CHR   Variant_ID nominal_P_value
#>   <chr> <chr>                <dbl>
#> 1 1     rs1            0.000000001
#> 2 1     rs2            0.05       
```
