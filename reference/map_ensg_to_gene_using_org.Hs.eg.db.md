# Map Ensembl IDs to Gene Symbols using `org.Hs.eg.db`

This function maps Ensembl IDs to their corresponding gene symbols using
the `org.Hs.eg.db` package.

## Usage

``` r
map_ensg_to_gene_using_org.Hs.eg.db(ensembl_ids, ensembl_col = NULL)
```

## Arguments

- ensembl_ids:

  A character vector of Ensembl IDs or a data frame containing Ensembl
  IDs.

- ensembl_col:

  The column name of Ensembl IDs if `ensembl_ids` is a data frame.

## Value

A data frame with two columns: the input Ensembl IDs and the
corresponding gene symbols.

## Examples

``` r
if (FALSE) { # \dontrun{
ensembl_ids <- c("ENSG00000141510.1", "ENSG00000012048", "ENSG00000146648") # Example with Ensembl ID vector
map_ensg_to_gene_using_org.Hs.eg.db(ensembl_ids = ensembl_ids)

ensembl_ids_df <- data.frame(EnsemblID = ensembl_ids, OtherInformation = c(1,2,3)) # Example with data frame input
map_ensg_to_gene_using_org.Hs.eg.db(ensembl_ids = ensembl_ids_df, ensembl_col = "EnsemblID")
} # }
```
