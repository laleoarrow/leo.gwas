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
if (FALSE) { # \dontrun{
# Suppose you have 25 elements and want to batch them in groups of 6
nums <- 1:25
it <- leo_iterator(nums, 6)
while (TRUE) {
  batch <- it()
  if (is.null(batch)) break

  # add your parallel process steps
  print(batch)
}
} # }
```
