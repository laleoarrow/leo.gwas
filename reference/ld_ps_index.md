# Loci_plot: Calculate the LD-matrix (LD r2) for the index SNP

Loci_plot: Calculate the LD-matrix (LD r2) for the index SNP

## Usage

``` r
ld_ps_index(
  gwas,
  index = "rs999",
  win = 1000,
  ld_calculation = T,
  bfile = "/Users/leoarrow/project/ref/1kg.v3/EAS",
  plink_bin = plinkbinr::get_plink_exe()
)
```

## Arguments

- gwas:

  gwas summary data that needs to select loci and calculate r2

- index:

  index snp

- win:

  window size to locally calculate the r2; set it larger than that you
  want to plot; \# default calculate 1MB

- ld_calculation:

  if calculate the LD locally, defaut T; use F if the index snp is a
  rare variant (MAF\<0.01)

- bfile:

  bfile

- plink_bin:

  plinkbinr::get_plink_exe()

## Value

loci data with calculated r2
