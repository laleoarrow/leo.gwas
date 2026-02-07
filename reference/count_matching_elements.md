# Count or Identify Matches of a Pattern in a Vector

This function counts the number of elements in a vector that contain a
given pattern, or returns a logical vector indicating which elements
match the pattern.

## Usage

``` r
count_matching_elements(vector, pattern, return = "count")
```

## Arguments

- vector:

  A character vector to be searched.

- pattern:

  The pattern to match (regular expression).

- return:

  A string specifying the return type: "count" for number of matches, or
  "logi" for a logical vector.

## Value

An integer representing the number of matches if return is "count", or a
logical vector indicating which elements match the pattern if return is
"logi".

## Examples

``` r
vec <- c("abc", "def", "xyz", "abcd")
count_matching_elements(vec, "abc", return = "count")  # Returns 2
#> [1] 2
count_matching_elements(vec, "abc", return = "logi")   # Returns c(TRUE, FALSE, FALSE, TRUE)
#> [1]  TRUE FALSE FALSE  TRUE
```
