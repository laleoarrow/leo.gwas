# Across a df to count TRUE and FALSE

This function summarizes both TRUE and FALSE values in each column of a
data frame.

## Usage

``` r
across_df_TF(df, type = "T")
```

## Arguments

- df:

  a data frame

- type:

  "T" for TRUE counts (default), "F" for FALSE counts

## Value

a data frame with the number of TRUE or FALSE values in each column

## Examples

``` r
df <- data.frame(a = c(TRUE, FALSE, TRUE, TRUE), b = c(FALSE, TRUE, TRUE, TRUE))
across_df_TF(df) # Count TRUE (default)
#>   a b
#> 1 3 3
across_df_TF(df, "F") # Count FALSE
#>   a b
#> 1 1 1
```
