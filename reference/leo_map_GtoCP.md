# Map Gene Symbols to Genomic Positions

TODO: Merge all gene annotation function into one simple command. This
function maps gene symbols to their genomic positions (chromosome,
start, end, strand) using the specified method and genome assembly.

## Usage

``` r
leo_map_GtoCP(
  genes,
  gene_col = NULL,
  method = c("bioconductor", "gtf"),
  genome = c("hg19", "hg38"),
  ...
)
```

## Arguments

- genes:

  A character vector of gene symbols or a data frame containing gene
  symbols.

- gene_col:

  The column name of gene symbols if `genes` is a data frame.

- method:

  The method to use: `"bioconductor"` or `"gtf"`.

- genome:

  The genome assembly to use: `"hg19"` or `"hg38"`.

- ...:

  Additional arguments to pass to the GTF method.

## Value

A data frame with columns: gene_symbol, chr, bp_start, bp_end, strand.

## Details

The function supports two methods:

- `"bioconductor"`: uses Bioconductor packages. See
  [`map_gene_to_chrbp_using_TxDb()`](https://laleoarrow.github.io/leo.gwas/reference/map_gene_to_chrbp_using_TxDb.md).

- `"gtf"`: uses a GTF file. See
  [`map_gene_to_chrbp_using_gtf()`](https://laleoarrow.github.io/leo.gwas/reference/map_gene_to_chrbp_using_gtf.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Using Bioconductor method with a character vector of gene symbols
leo_map_GtoCP(genes = c("TP53", "BRCA1", "EGFR"), method = "bioconductor", genome = "hg19")

# Using Bioconductor method with a data frame
leo_map_GtoCP(genes = data.frame(gene_name = c("TP53", "BRCA1", "EGFR"), value = c(1.2, 3.4, 5.6)),
              gene_col = "gene_name", method = "bioconductor", genome = "hg19")

# Using GTF method with a character vector of gene symbols
leo_map_GtoCP(genes = c("TP53", "BRCA1", "EGFR"), method = "gtf", genome = "hg38")

# Using GTF method with a data frame
leo_map_GtoCP(genes = data.frame(gene_name = c("TP53", "BRCA1", "EGFR"), value = c(1.2, 3.4, 5.6)),
              gene_col = "gene_name", method = "gtf", genome = "hg38", download_dir = "~/project/ref/gtf")
} # }
```
