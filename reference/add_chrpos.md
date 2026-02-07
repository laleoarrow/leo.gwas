# Convert rsID to CHR & BP

This function takes a dataset with a column containing rsIDs (SNP IDs)
and adds the corresponding chromosome (CHR) and position (POS)
information. It queries the SNPlocs.Hsapiens.dbSNP155.GRCh37 database
(or GRCh38 if specified) to retrieve the genomic positions. The function
returns a dataframe with the additional 'CHR' and 'POS' columns
appended.

## Usage

``` r
add_chrpos(dat, snp_col = "SNP", ref = "GRCh37")
```

## Arguments

- dat:

  A dataframe containing at least a column with SNP IDs (rsIDs).

- snp_col:

  A string indicating the column name containing SNP IDs (default is
  "SNP").

- ref:

  A string indicating the reference genome version. Default is "GRCh37",
  can also use "GRCh38".

## Value

A dataframe with additional 'CHR' and 'POS' columns.

## Examples

``` r
pacman::p_load(data.table, dplyr, BSgenome, SNPlocs.Hsapiens.dbSNP155.GRCh37)
#> Error in loadNamespace(x): there is no package called ‘pacman’
zuo_ref <- fread("/path/to/1KG-EAS-EAF.txt.gz") # Input dataset with rsID (SNP column)
#> Error in fread("/path/to/1KG-EAS-EAF.txt.gz"): could not find function "fread"
result <- add_chrpos(zuo_ref, snp_col = "SNP", ref = "GRCh37")
#> Error: object 'zuo_ref' not found
```
