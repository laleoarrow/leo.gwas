# Map Genes to Their TSS Positions

This function maps gene symbols to their transcription start sites (TSS)
positions using hg19 or hg38 genome assembly.

## Usage

``` r
map_gene_to_tss_using_gtf(
  genes,
  gene_col = NULL,
  genome = c("hg19", "hg38"),
  gtf_file = NULL,
  download_dir = "~/Project/ref/gtf",
  ...
)
```

## Arguments

- genes:

  A character vector of gene symbols or a data frame containing gene
  symbols.

- gene_col:

  The column name of gene symbols if `genes` is a data frame.

- genome:

  The genome assembly to use: `"hg19"` or `"hg38"`.

- gtf_file:

  The path to a GTF file. If `NULL`, the function will download the
  appropriate GTF file.

- download_dir:

  The path where you want to store the downloaded gtf file.

- ...:

  Pass to map_gene_to_chrbp_using_gtf. See
  [`map_gene_to_chrbp_using_gtf()`](https://laleoarrow.github.io/leo.gwas/reference/map_gene_to_chrbp_using_gtf.md)

## Value

A data frame with gene symbols and their TSS positions

## Examples

``` r
if (FALSE) { # \dontrun{
genes <- c("TP53", "BRCA1", "EGFR")
map_gene_to_tss_using_gtf(genes = genes, genome = "hg38")

genes_df <- data.frame(GeneName = genes, OtherInfo = c(1,2,3))
map_gene_to_tss_using_gtf(genes = genes_df, gene_col = "GeneName", genome = "hg19")
} # }
```
