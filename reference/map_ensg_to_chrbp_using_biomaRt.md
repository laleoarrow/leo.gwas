# Map Ensembl Gene IDs to Genomic Positions using biomaRt

This function uses biomaRt to retrieve genomic positions (chr, bp_start,
bp_end, and strand) from the Ensembl database based on Ensembl gene IDs.

## Usage

``` r
map_ensg_to_chrbp_using_biomaRt(
  ensembl_ids,
  ensembl_col = NULL,
  genome = c("hg19", "hg38"),
  verbose = F
)
```

## Arguments

- ensembl_ids:

  A character vector of Ensembl gene IDs, or a data frame containing
  this information.

- ensembl_col:

  If `ensembl_ids` is a data frame, specify the column name containing
  the Ensembl gene IDs.

- genome:

  The genome version to use: "hg19" or "hg38".

- verbose:

  Logical indicating whether to print the unmapped information.

## Value

A data frame containing genomic position information, including
ensembl_gene_id, chr, bp_start, bp_end, strand.

## Examples

``` r
if (FALSE) { # \dontrun{
# Query using Ensembl gene IDs
ensembl_ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648")
map_ensg_to_chrbp_using_biomaRt(ensembl_ids = ensembl_ids, genome = "hg19")

# Use a data frame as input
gene_df <- data.frame(ensembl_id = ensembl_ids, value = c(1.2, 3.4, 5.6))
map_ensg_to_chrbp_using_biomaRt(ensembl_ids = gene_df, ensembl_col = "ensembl_id", genome = "hg19")
} # }
```
