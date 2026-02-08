# Modified MR Scatter Plot

The `mr_scatter_plot` in the `TwoSampleMR` package could be better. This
is a modified version of `mr_scatter_plot` in the `TwoSampleMR` package
When the slope of two method is really close, the scatter plot may not
properly plot them! Change the line type here mannually herein then.

## Usage

``` r
mr_scatter_plot_modified(mr_results, dat)
```

## Arguments

- mr_results:

  same as TwoSampleMR::mr_scatter_plot

- dat:

  same as TwoSampleMR::mr_scatter_plot

## Examples

``` r
if (FALSE) { # \dontrun{
p1 <- mr_scatter_plot_modified(mr_results = res_pair, dat = dat_h_pair)
print(p1[[1]])
} # }
```
