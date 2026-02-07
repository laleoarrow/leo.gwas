# Loci_plot: save_regional_plot

Loci_plot: save_regional_plot

## Usage

``` r
save_regional_plot(
  path,
  loc,
  gene,
  save = T,
  title = expression(paste(italic("CLPSL1"), " (T1D)")),
  labels = c("index"),
  filter_gene_biotype = c("protein_coding"),
  border = F,
  width = 7.5,
  height = 5.5
)
```

## Arguments

- path:

  path to store the plot; make sure the path is exist

- loc:

  output from locuszoomr_loc

- gene:

  gene

- save:

  if T, will save plot to path; if F, return the plot only

- labels:

  labels; in case you need to indicate the index SNP and other SNP

- border:

  border for gene track

- width:

  width

- height:

  height
