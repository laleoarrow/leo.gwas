# Add rsID based on local reference file (recommended)

This function takes a data frame with at least chromosome, position, and
allele columns, and matches rsID based on a local reference file,
considering reversed/complement alleles. It allows for flexibility in
the names of the allele columns (e.g., 'REF'/'ALT' or 'A1'/'A2'). The
function returns the original data frame with rsID added in the first
column.

## Usage

``` r
add_rsid_using_ref(
  dat,
  local_ref,
  chr_col = "CHR",
  pos_col = "BP",
  a1_col = "A1",
  a2_col = "A2"
)
```

## Arguments

- dat:

  A data frame containing at least chromosome, position, and allele
  columns.

- local_ref:

  A data frame containing at least 'ID' and 'SNP' columns.

- chr_col:

  Name of the chromosome column in 'dat'. Default is 'CHR'.

- pos_col:

  Name of the position column in 'dat'. Default is 'BP'.

- a1_col:

  Name of the first allele column in 'dat'. Default is 'A1' (e.g., 'A1'
  or 'ALT').

- a2_col:

  Name of the second allele column in 'dat'. Default is 'A2' (e.g., 'A2'
  or 'REF').

## Value

The original data frame with an added 'rsID' column as the first column.

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
# Example data with REF and ALT
dat <- data.frame(
  CHR = c(7, 12, 4),
  BP = c(6013153, 126890980, 40088896),
  REF = c("G", "A", "T"),
  ALT = c("A", "G", "A")
)
# Reference data
local_ref <- data.frame(
  ID = c("7:6013153:A:G", "12:126890980:G:A", "4:40088896:A:T"),
  SNP = c("rs10000", "rs1000000", "rs10000000")
)
result <- add_rsid_using_ref(dat, local_ref, a1_col = "ALT", a2_col = "REF")
print(result)

# Example data with A1 and A2
dat2 <- data.frame(
  CHR = c(7, 12, 4),
  POS = c(6013153, 126890980, 40088896),
  A1 = c("A", "G", "A"),
  A2 = c("G", "A", "T")
)
result2 <- add_rsid_using_ref(dat2, local_ref, pos_col = "POS")
print(result2)
} # }
```
