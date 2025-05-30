% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smr_besd_tools.R
\name{leo_smr_extract_sig_res}
\alias{leo_smr_extract_sig_res}
\title{Extract significant results from \code{.all} files}
\usage{
leo_smr_extract_sig_res(dir, pass_type = c("FDR", "Bonferroni"), out_dir = "")
}
\arguments{
\item{dir}{Character. The main \code{combine_1outcome} folder from \func{combine_smr_res_1outcome}.}

\item{pass_type}{Character. Which significance criterion to use:
one of \code{c("FDR", "Bonferroni", "both")}.
- "FDR": keep rows where \code{Pass_FDR == "Pass"}
- "Bonferroni": keep rows where \code{Pass_Bonferroni == "Pass"}
- "both": keep rows where either FDR or Bonferroni is "Pass"}

\item{out_dir}{Character. Where to write significant results.
Default \code{dir/combine_1outcome/sig}.}
}
\value{
NULL. Writes \code{.sig.all} files to disk.
}
\description{
This function scans the \code{combine_1outcome} folder for \code{.all} files,
filters rows where \code{Pass_FDR == "Pass"} or \code{Pass_Bonferroni == "Pass"}
(depending on \code{pass_type} argument), and writes the significant subset
to a new file (e.g., \code{smr_res_YYYYMMDD_iri3.sig.all}).
}
\examples{
\dontrun{
# For the directory "/Users/leoarrow/project/iridocyclitis/output/smr-t2d",
# after running combine_smr_res_1outcome, we have .all files in
# "smr-t2d/combine_1outcome".

leo_smr_extract_sig_res(
  dir       = "/Users/leoarrow/project/iridocyclitis/output/smr-t2d",
  pass_type = "FDR"
)
}
}
