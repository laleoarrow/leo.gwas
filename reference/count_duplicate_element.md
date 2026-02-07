# Count or Identify Duplicates in a Vector

This function counts the number of duplicate elements in a vector, or
returns a logical vector indicating which elements are duplicates.

## Usage

``` r
count_duplicate_element(vector, return = "count")
```

## Arguments

- vector:

  A vector to be checked for duplicates.

- return:

  A string specifying the return type: "count" for number of duplicates,
  or "logi" for a logical vector.

## Value

An integer representing the number of duplicates if return is "count",
or a logical vector indicating which elements are duplicates if return
is "logi".

## Examples

``` r
vec <- c("a", "b", "c", "a", "b", "d")
count_duplicate_element(vec, return = "count")  # Returns 4
#> [1] 4
count_duplicate_element(vec, return = "logi")   # Returns c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
#> [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE
```
