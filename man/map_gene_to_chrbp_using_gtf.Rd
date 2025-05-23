% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0.annotation.R
\name{map_gene_to_chrbp_using_gtf}
\alias{map_gene_to_chrbp_using_gtf}
\title{Map Gene Symbols Using GTF File}
\usage{
map_gene_to_chrbp_using_gtf(
  genes,
  gene_col = NULL,
  genome = c("hg19", "hg38"),
  gtf_file = NULL,
  download_dir = "~/project/ref/gtf"
)
}
\arguments{
\item{genes}{A character vector of gene symbols or a data frame containing gene symbols.}

\item{gene_col}{The column name of gene symbols if \code{genes} is a data frame.}

\item{genome}{The genome assembly to use: \code{"hg19"} or \code{"hg38"}.}

\item{gtf_file}{The path to a GTF file. If \code{NULL}, the function will download the appropriate GTF file.}

\item{download_dir}{The path where you wanna store the downloaded gtf file}
}
\value{
A data frame with mapping results.
}
\description{
Map Gene Symbols Using GTF File
}
\examples{
\dontrun{
gene_symbols <- c("TP53", "BRCA1", "EGFR")
map_gene_to_chrbp_using_gtf(genes = gene_symbols, genome = "hg38")

gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInformation = c(1,2,3))
map_gene_to_chrbp_using_gtf(genes = gene_symbols_df, gene_col = "GeneName" , genome = "hg19")
}
}
