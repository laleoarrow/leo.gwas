# Map Gene Symbols to biotype & description via **annotables**

使用离线数据框
[`annotables::grch38`](https://rdrr.io/pkg/annotables/man/grch38.html)
与
[`annotables::grch37`](https://rdrr.io/pkg/annotables/man/grch37.html)
为基因符号添加 *biotype*、*description*，并给出来源列
`infer_version`（"grch38" / "grch37" / "unmapped"）。

## Usage

``` r
map_gene_class_using_annotables(genes, gene_col = "Gene", quiet = FALSE)
```

## Arguments

- genes:

  Character vector of gene symbols，或包含基因符号的数据框/ tibble

- gene_col:

  列名（当 `genes` 为表格时），默认 `"Gene"`

- quiet:

  逻辑值；`TRUE` 时不显示进度信息

## Value

与输入同结构的数据，附加 `biotype`、`description`、`infer_version`

## Details

逻辑：

1.  先左连接 **grch38**；命中则 `infer_version = "grch38"`

2.  未命中者再补 **grch37**；命中则 `infer_version = "grch37"`

3.  仍未命中者三列均填占位，如上所述

## Examples

``` r
if (FALSE) { # \dontrun{
map_gene_class_using_annotables(c("TP53","BRCA1","C14orf37"))

df <- tibble::tibble(Gene = c("TP53","XIST","C14orf37"))
map_gene_class_using_annotables(df, gene_col = "Gene")
} # }
```
