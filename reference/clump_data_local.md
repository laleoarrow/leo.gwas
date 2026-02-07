# Perform LD Clumping Locally or via Reference Panel

This function performs LD clumping using either a local PLINK binary
file (`bfile`) or a 1000 Genomes super population panel. Requires the
`ieugwasr` and `plinkbinr` packages.

## Usage

``` r
clump_data_local(
  dat,
  pop = NULL,
  bfile = NULL,
  clump_kb = 10000,
  clump_r2 = 0.001,
  plink_bin = plinkbinr::get_plink_exe()
)
```

## Arguments

- dat:

  Data frame with columns `SNP`, `pval.exposure`; optional
  `id.exposure`.

- pop:

  1000G super-pop for online clumping (AFR/AMR/EAS/EUR/SAS). Used when
  `bfile` is NULL.

- bfile:

  PLINK reference panel prefix for local clumping (without
  .bed/.bim/.fam). If set, `pop` is ignored.

- clump_kb:

  Window size in kb. Default 10000.

- clump_r2:

  r^2 threshold. Default 0.001.

- plink_bin:

  Path to PLINK executable (auto-detect via plinkbinr if NULL and
  `bfile` is set).

## Examples

``` r
# Reference:
# - https://github.com/MRCIEU/TwoSampleMR/issues/173
# - https://blog.csdn.net/xiaozheng1213/article/details/126269969
library(ieugwasr)
#> OpenGWAS updates:
#>   Date: 2024-05-17
#>   [>] OpenGWAS is growing!
#>   [>] Please take 2 minutes to give us feedback -
#>   [>] It will help directly shape our emerging roadmap
#>   [>] https://forms.office.com/e/eSr7EFAfCG
library(plinkbinr) # devtools::install_github("explodecomputer/plinkbinr")
#> Error in library(plinkbinr): there is no package called ‘plinkbinr’
plinkbinr::get_plink_exe()
#> Error in loadNamespace(x): there is no package called ‘plinkbinr’

# Note: after using this, please check `leo_clump` column to see if they are all TRUE
# If it's contains F, it means no SNPs remained after clumping or something bad happened
```
