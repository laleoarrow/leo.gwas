# Map Gene Symbols Using GTF File

Map Gene Symbols Using GTF File

## Usage

``` r
map_gene_to_chrbp_using_gtf(
  genes,
  gene_col = NULL,
  genome = c("hg19", "hg38"),
  gtf_file = NULL,
  download_dir = "~/project/ref/gtf"
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

  The path where you wanna store the downloaded gtf file

## Value

A data frame with mapping results.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_symbols <- c("TP53", "BRCA1", "EGFR")
map_gene_to_chrbp_using_gtf(genes = gene_symbols, genome = "hg38")

gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInformation = c(1,2,3))
map_gene_to_chrbp_using_gtf(genes = gene_symbols_df, gene_col = "GeneName" , genome = "hg19")
} # }
```
