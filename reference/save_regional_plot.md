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

  Path to the directory containing summary statistics.

- loc:

  Locus string (e.g., "1:1000-2000").

- gene:

  Gene name to highlight.

- save:

  Logical. Whether to save the plot (default: TRUE).

- title:

  Title of the plot.

- labels:

  Labels to indicate index SNP and other SNPs.

- filter_gene_biotype:

  Character vector of gene biotypes to filter.

- border:

  Logical. Whether to add a border for gene track.

- width:

  Plot width.

- height:

  Plot height.
