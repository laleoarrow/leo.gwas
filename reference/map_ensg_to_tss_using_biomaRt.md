# Map Ensembl Gene IDs to TSS Using biomaRt

This function maps Ensembl gene IDs to their transcription start sites
(TSS). (Note that this is inferred using the gene start and end
information based on strand) It uses biomaRt for the specified genome
assembly ("hg19" or "hg38").

## Usage

``` r
map_ensg_to_tss_using_biomaRt(
  ensembl_ids,
  ensembl_col = NULL,
  genome = c("hg19", "hg38"),
  ...
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

- ...:

  pass to `map_ensg_to_chrbp_using_biomaRt`. See
  [`map_ensg_to_chrbp_using_biomaRt()`](https://laleoarrow.github.io/leo.gwas/reference/map_ensg_to_chrbp_using_biomaRt.md)

## Value

A data frame with Ensembl gene IDs and their TSS positions.

## Examples

``` r
if (FALSE) { # \dontrun{
ensembl_ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648")
map_ensg_to_tss_using_biomaRt(ensembl_ids = ensembl_ids, genome = "hg19")

ensembl_ids_df <- data.frame(EnsemblID = ensembl_ids, OtherInfo = c(1, 2, 3))
map_ensg_to_tss_using_biomaRt(ensembl_ids = ensembl_ids_df,
                              ensembl_col = "EnsemblID",
                              genome = "hg38")
} # }
```
