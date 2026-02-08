# Map Gene Symbols to biotype & description via **annotables**

Use local dataframes
[`annotables::grch38`](https://rdrr.io/pkg/annotables/man/grch38.html)
and
[`annotables::grch37`](https://rdrr.io/pkg/annotables/man/grch37.html)
to add *biotype* and *description* to gene symbols, providing a source
column `infer_version` ("grch38" / "grch37" / "unmapped").

## Usage

``` r
map_gene_class_using_annotables(genes, gene_col = "Gene", quiet = FALSE)
```

## Arguments

- genes:

  Character vector of gene symbols, or data frame/tibble containing gene
  symbols

- gene_col:

  Column name (when `genes` is a table), default `"Gene"`

- quiet:

  Logical; if `TRUE`, suppress progress messages

## Value

Data frame with same structure as input, plus `biotype`, `description`,
`infer_version`

## Details

Logic:

1.  Left join with **grch38**; if matched, `infer_version = "grch38"`

2.  If not matched, fill with **grch37**; if matched,
    `infer_version = "grch37"`

3.  If still not matched, fill columns with placeholders as described
    above

## Examples

``` r
if (FALSE) { # \dontrun{
map_gene_class_using_annotables(c("TP53","BRCA1","C14orf37"))

df <- tibble::tibble(Gene = c("TP53","XIST","C14orf37"))
map_gene_class_using_annotables(df, gene_col = "Gene")
} # }
```
