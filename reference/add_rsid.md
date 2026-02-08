# Convert CHR:BP to rsID (not recommended)

This function takes a dataframe that must include columns 'CHR' and
'BP', and it appends the corresponding rsID by querying the
SNPlocs.Hsapiens.dbSNP155.GRCh37 database. The function returns a
dataframe with the rsIDs included.

## Usage

``` r
add_rsid(dat, ref = "GRCh37")
```

## Arguments

- dat:

  - "GRCh38" for GRCh38 Ref panel

- ref:

  str indicating ref version.

  - "GRCh37" for GRCh37 Ref panel

## Value

A dataframe with an additional 'RefSNP_id' column which contains the
rsIDs.

## Examples

``` r
pacman::p_load(data.table, BSgenome, leo.gwas)
#> Error in loadNamespace(x): there is no package called 'pacman'
library("SNPlocs.Hsapiens.dbSNP155.GRCh37") # for GRCh37
#> Error in library("SNPlocs.Hsapiens.dbSNP155.GRCh37"): there is no package called 'SNPlocs.Hsapiens.dbSNP155.GRCh37'
library("SNPlocs.Hsapiens.dbSNP155.GRCh38") # for GRCh38
#> Error in library("SNPlocs.Hsapiens.dbSNP155.GRCh38"): there is no package called 'SNPlocs.Hsapiens.dbSNP155.GRCh38'
df <- data.frame(
  CHR = c(1, 1),
  BP = c(15211, 15820)
)
result <- add_rsid(df); result
#> Error in loadNamespace(x): there is no package called 'SNPlocs.Hsapiens.dbSNP155.GRCh37'
```
