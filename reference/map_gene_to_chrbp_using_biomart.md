# Map Gene Symbols to Genomic Positions Using biomaRt

This function queries gene symbols for their genomic positions
(chromosome, start, end, strand) using Ensembl's biomaRt for the
specified genome assembly ("hg19" or "hg38").

## Usage

``` r
map_gene_to_chrbp_using_biomaRt(
  genes,
  gene_col = NULL,
  genome = c("hg19", "hg38")
)
```

## Arguments

- genes:

  A character vector of gene symbols to query, or a data frame
  containing gene symbols.

- gene_col:

  The column name of gene symbols if `genes` is a data frame.

- genome:

  The genome assembly to use: "hg19" or "hg38".

## Value

A data frame with the original data and new columns: chr, bp_start,
bp_end, strand, gene_symbol.

## Details

If the input is a data frame, the function retains all existing columns
and adds new columns with the mapping results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Query location of TP53, BRCA1, and EGFR genes
gene_symbols <- c("TP53", "BRCA1", "EGFR")
map_gene_to_chrbp_using_biomaRt(genes = gene_symbols, genome = "hg19")

# Using a data frame with gene symbols
gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInfo = c(1, 2, 3))
map_gene_to_chrbp_using_biomaRt(genes = gene_symbols_df, gene_col = "GeneName", genome = "hg19")
} # }
```
