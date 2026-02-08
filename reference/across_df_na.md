# Across a df to count na

This function summarize NA values in each column of a data frame.

## Usage

``` r
across_df_na(df)
```

## Arguments

- df:

  a data frame

## Value

a data frame with the number of NA values in each column

## Examples

``` r
df <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 2, 3, 4))
across_df_na(df)
#>   a b
#> 1 1 1
```
