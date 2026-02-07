# Map Gene Symbols to Ensembl IDs

This function provides robust gene symbol to Ensembl ID mapping through:

1.  Local `org.Hs.eg.db` annotations (default)

2.  Ensembl BioMart web service (requires internet)

## Usage

``` r
map_gene_to_ensembl(
  genes,
  gene_col = NULL,
  method = c("org.Hs.eg.db", "biomart"),
  genome = c("hg19", "hg38"),
  type = c("first", "combine"),
  sep = "/",
  batch_size = 100
)
```

## Arguments

- genes:

  Input containing gene symbols. Can be:

  - Character vector of gene symbols

  - Data frame containing gene symbol column

- gene_col:

  Column name containing gene symbols when `genes` is data frame

- method:

  Mapping methodology:

  - "org.Hs.eg.db": Local Bioconductor annotations (default)

  - "biomart": Ensembl BioMart service

- genome:

  Genome assembly version (BioMart only):

  - "hg19": GRCh37 (default)

  - "hg38": GRCh38

- type:

  Multi-mapping handling:

  - "first": Return first valid ID (default)

  - "combine": Concatenate multiple IDs

- sep:

  Separator for combined IDs (default: "/")

- batch_size:

  BioMart query batch size (default: 100)

## Value

Data frame with original data + ensembl_id column

## Examples

``` r
if (FALSE) { # \dontrun{
# Local annotation method
genes <- c("TP53", "BRCA1", "VEGFA")
result_local <- map_gene_to_ensembl(genes)

# BioMart with custom parameters
gene_df <- data.frame(
  my_symbol = c("TP53", "BRCA1", "NONEXISTENT"),
  values = rnorm(3)
)
result_biomart <- map_gene_to_ensembl(
  gene_df,
  gene_col = "my_symbol",
  method = "biomart",
  genome = "hg19",
  batch_size = 50
)
} # }
```
