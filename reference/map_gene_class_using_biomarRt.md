# Map Gene Symbols to biotype & description via **biomaRt** (GRCh38 -\> GRCh37 fallback)

First connect to **Ensembl GRCh38** (www.ensembl.org) to query
*gene_biotype* and *description* in batch; for unmapped genes,
automatically fallback to **GRCh37** (grch37.ensembl.org). Finally
append three columns for each gene:

- **biotype**

- **description**

- **infer_version** – "GRCh38" / "GRCh37" / "unmapped"

## Usage

``` r
map_gene_class_using_biomarRt(genes, gene_col = "Gene", quiet = FALSE)
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

## Examples

``` r
if (FALSE) { # \dontrun{
map_gene_class_using_biomarRt(c("TP53", "IGHV1-69", "TRAC", "MIR21"))

df <- tibble::tibble(Gene = c("TP53","XIST","SNORD14A"))
map_gene_class_using_biomarRt(df, gene_col = "Gene")
} # }
```
