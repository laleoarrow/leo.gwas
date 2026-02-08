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
p1 <- mr_scatter_plot_modified(mr_results = res_pair, dat = dat_h_pair)
#> Error: object 'dat_h_pair' not found
print(p1[[1]])
#> Error: object 'p1' not found
library(ggplot2);library(ggsci);library(TwoSampleMR);library(dplyr)
#> TwoSampleMR version 0.6.30 
#>   [>] New authentication requirements: https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication.
#>   [>] Major upgrades to our servers completed to improve service and stability.
#>   [>] We need your help to shape our emerging roadmap!
#>       Please take 2 minutes to give us feedback -
#>       https://forms.office.com/e/eSr7EFAfCG
#> 
#> Attaching package: 'TwoSampleMR'
#> The following object is masked from 'package:ieugwasr':
#> 
#>     ld_matrix
```
