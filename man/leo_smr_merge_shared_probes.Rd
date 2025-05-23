% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smr_besd_tools.R
\name{leo_smr_merge_shared_probes}
\alias{leo_smr_merge_shared_probes}
\title{Merge multiple SMR files and keep only shared probes for 2 outcomes}
\usage{
leo_smr_merge_shared_probes(
  dir = NULL,
  file_paths = NULL,
  pattern = "\\\\.sig\\\\.all$",
  out_file = ""
)
}
\arguments{
\item{dir}{Character. A directory containing SMR files (e.g. ".sig.all" files).
If provided, the function will read all matching files in this directory.}

\item{file_paths}{Character vector. A set of absolute paths to SMR files. (Priority over \code{dir})
If this is non-NULL, \code{dir} is ignored.}

\item{pattern}{Character. Regex pattern for searching files in \code{dir}. Default is "\\.sig\\.all$".}

\item{out_file}{Character. If not empty, write the merged result to this path. Otherwise only return in R.}
}
\value{
A tibble containing shared probes across all specified files. If \code{out_file} is provided, the result is written to disk and the function returns \code{NULL}.
}
\description{
This function finds intersected \code{probeID} across all provided files
and merges them into a single data frame using \code{inner_join}.
}
\examples{
\dontrun{
# 1) Merge all .sig.all files in a directory:
merged_df <- leo_smr_merge_shared_probes(
  dir = "/path/to/smr-t2d/combine_1outcome/sig"
)

# 2) Merge specific files by absolute paths:
files_vec <- c(
  "/path/to/smr/combine_1outcome/smr_res_20241230_iri3.all",
  "/path/to/smr/combine_1outcome/smr_res_20241230_t1d1.all"
)
merged_df <- leo_smr_merge_shared_probes(
  file_paths = files_vec,
  out_file = "shared_probes_merged.tsv"
)
}
}
