# Map Gene Symbols Using Bioconductor Packages

Map Gene Symbols Using Bioconductor Packages

## Usage

``` r
map_gene_to_chrbp_using_TxDb(
  genes,
  gene_col = NULL,
  genome = c("hg19", "hg38")
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

  - For hg19, it needs `"TxDb.Hsapiens.UCSC.hg19.knownGene"`
    Bioconductor package

  - For hg38, it needs `"TxDb.Hsapiens.UCSC.hg38.knownGene"`
    Bioconductor package

## Value

A data frame with mapped results.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_symbols <- c("TP53", "BRCA1", "EGFR") # Example with gene symbol vector
map_gene_to_chrbp_using_TxDb(genes = gene_symbols, genome = "hg19")

gene_symbols_df <- data.frame(GeneName = gene_symbols,
                            OtherInformation = c(1,2,3)) # Example with data frame input
map_gene_to_chrbp_using_TxDb(genes = gene_symbols_df,
                             gene_col = "GeneName",
                             genome = "hg19")
} # }
```
