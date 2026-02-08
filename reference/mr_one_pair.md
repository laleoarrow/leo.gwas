# One-Click Perform 2SMR

`perform_mr_for_one_pair` perform a MR for one pair of exp and out

## Usage

``` r
mr_one_pair(
  dat_h,
  exp = "",
  out = "",
  save_plot = T,
  res_dir = "./output/tsmr",
  fig_dir = "./figure/tsmr"
)
```

## Arguments

- dat_h:

  harmonized data

- exp:

  a str indicating the exposure in dat_h

- out:

  a str indicating the outcome in dat_h

- save_plot:

  only save the plot if T, defaut T.

- res_dir:

  dir path where the result fo MR analysis stored

- fig_dir:

  dir path where the figure of MR analysis stored

## Value

list(res_pair=res_pair, res_pair_presso=res_pair_presso)

## Examples

``` r
if (FALSE) { # \dontrun{
# dat_h is harmonized TwoSampleMR data with columns:
# exposure, outcome, mr_keep, beta.exposure, beta.outcome, etc.
out <- mr_one_pair(
  dat_h = dat_h,
  exp = "BMI",
  out = "T1D",
  save_plot = FALSE,
  res_dir = "./output/tsmr",
  fig_dir = "./figure/tsmr"
)
names(out)
} # }
```
