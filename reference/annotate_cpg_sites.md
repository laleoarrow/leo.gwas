# Annotate CpG Sites with Gene Information

This function annotates a vector of CpG site probe IDs by retrieving
corresponding gene names from the Illumina 450k annotation package. It
returns a data frame with the original CpG sites and their associated
gene names.

## Usage

``` r
annotate_cpg_sites(
  cpg_vector,
  annotation_package = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)
```

## Arguments

- cpg_vector:

  A character vector of CpG site probe IDs to be annotated.

- annotation_package:

  A character string specifying the Illumina annotation package to use.
  It can be one of the following:

  - "IlluminaHumanMethylation450kanno.ilmn12.hg19" (default)

  - "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"

## Value

A data frame with two columns:

- `CpG_Site`: The original CpG site probe IDs.

- `Gene`: The associated gene names. `NA` if no gene annotation is
  found.

## Examples

``` r
if (FALSE) { # \dontrun{
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# Example CpG probe IDs
cpg_ids <- c("cg00000029", "cg00000108", "cg00000109")
# Annotate CpG sites
annotate_cpg_sites(cpg_ids, "IlluminaHumanMethylation450kanno.ilmn12.hg19")
} # }
```
