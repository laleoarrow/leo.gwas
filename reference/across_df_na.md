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
summarize_na(df)
#> Error in summarize_na(df): could not find function "summarize_na"
```
