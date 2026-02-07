# Helper functions for `leo.gwas_qc`

Internal utilities used in the GWAS QC pipeline.

## Usage

``` r
is_complementary(a1st, a2nd)

fetch_indel(df, type = "both")

fetch_non_indel(df)

duplicated_SNP_lines(
  df,
  type = "rm",
  dup_columns = c("SNP"),
  group_columns = dup_columns
)

slice1_SNP_lines(df, dup_columns = c("SNP"), group_columns = dup_columns)

fetch_same_direcrtion(df_x2, df_lg)

any_na(df)
```

## Functions

- `is_complementary(a1st, a2nd)`:

  Check if two alleles form an A/T or C/G pair.

- `fetch_indel(df, type)`:

  Filter indels by allele-string length.

- `fetch_non_indel(df)`:

  Keep SNPs with single-base alleles.

- `duplicated_SNP_lines(df, type, dup_columns, group_columns)`:

  Get/remove duplicated SNP rows.

- `slice1_SNP_lines(df, dup_columns, group_columns)`:

  Within duplicated groups, keep first row.

- `fetch_same_direcrtion(df_x2, df_lg)`:

  Keep same-direction effects between datasets.

- `any_na(df)`:

  Count NAs per column.

## Value

- `is_complementary`: logical vector.

- `fetch_indel`, `fetch_non_indel`, `slice1_SNP_lines`,
  `duplicated_SNP_lines("rm")`: data frame.

- `duplicated_SNP_lines("get")`: data frame with `count` column.

- `fetch_same_direcrtion`: filtered data frame (`df_x2` subset).

- `any_na`: named numeric vector of NA counts.

## Examples

``` r
# Small demo dataset
df <- data.frame(
  SNP = c("rs1","rs2","rs3","rs2"),
  A1  = c("A","AT","C","C"),
  A2  = c("T","A","G","G"),
  OR  = c(1.2, 0.9, 1.1, 1.1)
)

# Complementary alleles
is_complementary("A","T")  # TRUE
#> [1] TRUE
is_complementary("A","G")  # FALSE
#> [1] FALSE

# Indels and non-indels
fetch_indel(df, "both")
#> [1] SNP A1  A2  OR 
#> <0 rows> (or 0-length row.names)
fetch_non_indel(df)
#>   SNP A1 A2  OR
#> 1 rs1  A  T 1.2
#> 2 rs3  C  G 1.1
#> 3 rs2  C  G 1.1

# Duplicates by SNP
duplicated_SNP_lines(df, "get", dup_columns = "SNP")
#> Duplicated SNP line count table
#> 
#> 2 
#> 2 
#> # A tibble: 2 × 5
#>   SNP   A1    A2       OR count
#>   <chr> <chr> <chr> <dbl> <int>
#> 1 rs2   AT    A       0.9     2
#> 2 rs2   C     G       1.1     2
slice1_SNP_lines(df, dup_columns = "SNP")
#>   SNP A1 A2  OR
#> 1 rs1  A  T 1.2
#> 2 rs2 AT  A 0.9
#> 3 rs3  C  G 1.1

# Same-direction effects between two datasets
df2 <- transform(df, OR = c(1.1, 1.3, 0.8, 1.1))
fetch_same_direcrtion(df, df2)
#> The number with opposite effect size：2
#>   SNP A1 A2  OR
#> 1 rs1  A  T 1.2
#> 2 rs2  C  G 1.1

# Count NA by column
any_na(df)
#> SNP  A1  A2  OR 
#>   0   0   0   0 
```
