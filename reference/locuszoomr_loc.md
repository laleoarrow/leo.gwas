# Loci_plot: prepare the locus data for locuszoomr

Loci_plot: prepare the locus data for locuszoomr

## Usage

``` r
locuszoomr_loc(loci_data, gene, online_ld = F, index_snp, flank)
```

## Arguments

- loci_data:

  output from ld_ps_index

- gene:

  gene loci

- online_ld:

  whether to use online LD; default is F

- index_snp:

  indexed snp

- flank:

  flank size for the locus plot

## Value

prepared data which could be pass to save_regional_plot
