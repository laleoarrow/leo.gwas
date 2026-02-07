# Leo batch iterator builder

Sometimes you need to handle a large number of iterations, but
multi-core parallel computing can be tricky—memory usage may grow beyond
what is actually required. This function helps by batching the
iterations, so you only process a limited number of elements at a time,
preventing excessive RAM consumption.

## Usage

``` r
leo_iterator(elements, batch_size)
```

## Arguments

- elements:

  A vector or list to be iterated.

- batch_size:

  Number of elements per batch.

## Value

A function (no arguments). Repeated calls produce successive batches or
`NULL` if finished.

## Details

Each call to the returned function yields the next batch. Returns `NULL`
when no more elements remain.

## Examples

``` r
# Suppose you have 25 elements and want to batch them in groups of 6
nums <- 1:25
it <- leo_iterator(nums, 6)
#> ℹ [18:09:05] Creating batch iterator with 25 elements in total;  batch size = 6 ; total rounds = 5 
while (TRUE) {
  batch <- it()
  if (is.null(batch)) break

  # add your parallel process steps
  print(batch)
}
#> ✔ [18:09:05] Round 1 of 5 | elements from index 1 to 6 
#> [1] 1 2 3 4 5 6
#> ✔ [18:09:05] Round 2 of 5 | elements from index 7 to 12 
#> [1]  7  8  9 10 11 12
#> ✔ [18:09:05] Round 3 of 5 | elements from index 13 to 18 
#> [1] 13 14 15 16 17 18
#> ✔ [18:09:05] Round 4 of 5 | elements from index 19 to 24 
#> [1] 19 20 21 22 23 24
#> ✔ [18:09:05] Round 5 of 5 | elements from index 25 to 25 
#> [1] 25
```
