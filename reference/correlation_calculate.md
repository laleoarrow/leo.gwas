# Calculate Correlation between Two Vectors

This function calculates the Spearman (default) or Pearson correlation
coefficient and its associated p-value between two vectors. It
automatically handles missing values.

## Usage

``` r
correlation_calculate(vector_x, vector_y, method = "spearman", ...)
```

## Arguments

- vector_x:

  A numeric vector.

- vector_y:

  A numeric vector of the same length as `vector_x`.

- method:

  A character string specifying the correlation method ("spearman" or
  "pearson"). Defaults to "spearman".

- ...:

  Pass to [`cor.test`](https://rdrr.io/r/stats/cor.test.html).

## Value

A list with the correlation coefficient and p-value.

## See also

[`correlation_draw`](https://laleoarrow.github.io/leo.gwas/reference/correlation_draw.md)
for plotting correlation results.

## Examples

``` r
vector_x <- c(1, 2, 2, 4, 5)
vector_y <- c(5, 6, 7, 8, 7)
result <- correlation_calculate(vector_x, vector_y, method = "pearson")
```
