# Give precise p-value for chi-square test

Sometimes you get a p-value of 0 when you perform a chi-square test or
other analysis. This is because the p-value is so small that R rounds it
to 0. This function gives you a more precise p-value.

## Usage

``` r
chisq_p_value(chisq_value, df, digits = 4, prec = 100)
```

## Arguments

- chisq_value:

  numeric, chi-square value

- df:

  numeric, degree of freedom

- digits:

  numeric, digits for output to illustrate

- prec:

  numeric, precision for mpfr() function

## Value

A precise p-value in scientific format

## Examples

``` r
# install.packages("Rmpfr")
library(Rmpfr)
#> Error in library(Rmpfr): there is no package called 'Rmpfr'
chisq_value <- 2629; df <- 1
p_value <- chisq_p_value(chisq_value, df)
#> Error in mpfr(log_p_value, prec = prec): could not find function "mpfr"
print(p_value)
#> Error: object 'p_value' not found
```
