# Map Gene Symbols to biotype & description via **biomaRt** (GRCh38 → GRCh37 fallback)

先连接 **Ensembl GRCh38**（www.ensembl.org）批量查询 *gene_biotype* 与
*description*；未命中的基因再自动回退到
**GRCh37**（grch37.ensembl.org）。 最终为每个基因追加三列：

- **biotype**

- **description**

- **infer_version** — "GRCh38" / "GRCh37" / "unmapped"

## Usage

``` r
map_gene_class_using_biomarRt(genes, gene_col = "Gene", quiet = FALSE)
```

## Arguments

- genes:

  Character vector of gene symbols，或包含基因符号的数据框 / tibble

- gene_col:

  列名（当 `genes` 为表格时），默认 `"Gene"`

- quiet:

  逻辑值；`TRUE` 时不显示进度信息

## Value

输入同结构的数据，附加 `biotype`、`description`、`infer_version`

## Examples

``` r
if (FALSE) { # \dontrun{
map_gene_class_using_biomarRt(c("TP53", "IGHV1-69", "TRAC", "MIR21"))

df <- tibble::tibble(Gene = c("TP53","XIST","SNORD14A"))
map_gene_class_using_biomarRt(df, gene_col = "Gene")
} # }
```
