# Get Unique Identifier for Genetic Data

Get ID using CHR, BP, A2 (REF/Non-effect), A1 (ALT/Effect)

## Usage

``` r
get_id(x, count_A1_A2 = FALSE)
```

## Arguments

- x:

  A data.frame that must contain the columns CHR, BP, A2, and A1. Each
  column represents:

  - CHR: Chromosome (It can be any in c("chrom", "CHR", "Chromosome",
    "chromosome", "Chr"))

  - BP/POS: Base pair position (It can be any in c("pos", "POS",
    "position", "BP", "Position", "Bp"))

  - A2: Reference allele/non-effect allele (It can be any in c("A2",
    "Allele2", "allele2", "a2", "REF", "Ref", "ref", "Non-effect"))

  - A1: Alternative allele/effect allele (It can be any in c("A1",
    "Allele1", "allele1", "a1", "ALT", "Alt", "alt", "Effect"))

- count_A1_A2:

  If TRUE, will count the number of characters in A1 and A2

## Value

A data.frame with an additional 'ID' column (if count_A1_A2=TRUE,
containing unique identifiers and character counts of A1 and A2)

## Examples

``` r
df <- data.frame(chrom = c(1, 1, 2), pos = c(12345, 54321, 11111),
                 A1 = c("A", "T", "G"), A2 = c("G", "C", "A"))
get_id(df); get_id(df, count_A1_A2 = TRUE)
#>   chrom   pos A1 A2          ID
#> 1     1 12345  A  G 1:12345:G:A
#> 2     1 54321  T  C 1:54321:C:T
#> 3     2 11111  G  A 2:11111:A:G
#>   chrom   pos A1 A2          ID A1_n A2_n
#> 1     1 12345  A  G 1:12345:G:A    1    1
#> 2     1 54321  T  C 1:54321:C:T    1    1
#> 3     2 11111  G  A 2:11111:A:G    1    1
```
