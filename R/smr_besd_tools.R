#' This contains need-to-have functions for making BESD files for SMR and HEIDI



#' Filter Chromosomes Based on SNP P-value Threshold
#'
#' This function adds an index column to the input data frame and filters chromosomes based on whether any SNP within the chromosome crosses a specified threshold. If a chromosome has at least one SNP that meets the threshold, all SNPs on that chromosome are retained. Otherwise, all SNPs on that chromosome are removed.
#'
#' @param df A data frame containing SNP data.
#' @param chr_col Character string specifying the name of the chromosome column. Default is `"CHR"`.
#' @param snp_col Character string specifying the name of the SNP identifier column. Default is `"Variant_ID"`.
#' @param p_val_col Character string specifying the name of the p-value column. Default is `"nominal_P_value"`.
#' @param threshold Numeric value specifying the p-value threshold. Default is `5e-8`.
#'
#' @return A filtered data frame with an added `index` column.
#' @export
#' @examples
#' library(dplyr)
#' eqtl_data <- data.frame(
#' CHR = c("1", "1", "2", "2", "3"),
#' Variant_ID = c("rs1", "rs2", "rs3", "rs4", "rs5"),
#' nominal_P_value = c(1e-9, 0.05, 0.2, 1e-7, 0.3)
#' );eqtl_data
#' df_filtered <- filter_chr_basedonSNP_p(
#'   df = eqtl_data,
#'   chr_col = "CHR",
#'   snp_col = "Variant_ID",
#'   p_val_col = "nominal_P_value",
#'   threshold = 5e-8
#' );df_filtered
filter_chr_basedonSNP_p <- function(df,
                               chr_col = "CHR",
                               snp_col = "Variant_ID",
                               p_val_col = "nominal_P_value",
                               threshold = 1.57e-3) {
  chr_sym <- sym(chr_col)
  snp_sym <- sym(snp_col)
  pval_sym <- sym(p_val_col)

  df_filtered <- df %>%
    dplyr::group_by(!!chr_sym) %>%
    dplyr::filter(any((!!pval_sym) < threshold)) %>%
    dplyr::ungroup()
  return(df_filtered)
}

#' Filter Chromosomes Based on SNP P-value Threshold for `.qtltoolsnomi` files
#'
#' This function adds an index column to the input data frame and filters chromosomes based on whether any SNP within the chromosome crosses a specified threshold. If a chromosome has at least one SNP that meets the threshold, all SNPs on that chromosome are retained. Otherwise, all SNPs on that chromosome are removed.
#'
#' @param df A data frame containing SNP data.
#' @param chr_col Character string specifying the name of the chromosome column. Default is `"CHR"`.
#' @param snp_col Character string specifying the name of the SNP identifier column. Default is `"Variant_ID"`.
#' @param gene_col Character string specifying the name of the gene column. Default is `"Gene"`.
#' @param p_val_col Character string specifying the name of the p-value column. Default is `"nominal_P_value"`.
#' @param threshold Numeric value specifying the p-value threshold. Default is `5e-8`.
#'
#' @return A filtered data frame with an added `index` column.
#' @export
#' @examples
#' library(dplyr)
#' eqtl_data <- data.frame(
#' Gene = c("G1", "G1", "G2", "G2", "G3"),
#' CHR = c("1", "1", "2", "2", "3"),
#' Variant_ID = c("rs1", "rs2", "rs3", "rs4", "rs5"),
#' nominal_P_value = c(1e-9, 0.05, 0.2, 1e-7, 0.3)
#' );eqtl_data
#' df_filtered <- filter_chr_basedonSNP_p_qtltools(
#'   df = eqtl_data,
#'   chr_col = "CHR",
#'   snp_col = "Variant_ID",
#'   p_val_col = "nominal_P_value",
#'   threshold = 5e-8
#' );df_filtered
filter_chr_basedonSNP_p_qtltools <- function(df,
                                    chr_col = "CHR",
                                    snp_col = "Variant_ID",
                                    gene_col = "Gene",
                                    p_val_col = "nominal_P_value",
                                    threshold = 1.57e-3) {
  chr_sym <- sym(chr_col)
  snp_sym <- sym(snp_col)
  pval_sym <- sym(p_val_col)
  gene_sym <- sym(gene_col)

  df_filtered <- df %>%
    dplyr::group_by(!!gene_sym, !!chr_sym) %>%
    dplyr::filter(any((!!pval_sym) < threshold)) %>%
    dplyr::ungroup()
  return(df_filtered)
}


# Below is how you aggregate all SMR results for XWAS ----
# smr result should be all stored in a dir in the following format:
# - level 1: qtl type (e.g. mqtl, eqtl, ...)
# - level 2: qtl specific source (e.g. GTEx49, ...)

# 1. deal with single SMR results from one source ----
#' Combine SMR Results for All Chromosomes
#'
#' This function combines SMR (Summary-data-based Mendelian Randomization) result files across all chromosomes
#' for each unique exposure and outcome pair. The combined results are saved to a specified output directory.
#'
#' @param dir Character. The directory containing SMR result files. Files should follow the naming convention
#'            `exposure_chrX@outcome.smr`, where `X` represents the chromosome number.
#' @param out_dir Character. The output directory where combined SMR files will be saved.
#'                If not specified, defaults to a subdirectory named `chr_combined` within `dir`.
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Combine SMR results in the "data/smr_results" directory and save to default output directory
#' combine_smr_res_chr(dir = "data/smr_results")
#'
#' # Combine SMR results and specify a custom output directory
#' combine_smr_res_chr(dir = "data/smr_results", out_dir = "data/combined_results")
#' }
#'
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning cli_alert_danger
#' @importFrom dplyr filter pull
#' @importFrom data.table fread fwrite
#' @export
combine_smr_res_chr <- function(dir, out_dir="") {
  leo.gwas::leo_log("Combine SMR results for all chromosomes.")
  df_tmp <- dplyr::tibble(
    files=list.files(dir, full.names = F, pattern = "smr$"),
    exposures=strsplit(files, "@") %>% sapply(function(x) x[1]) %>% sub("_chr[0-9]+", "", .),
    outcomes=strsplit(files, "@") %>% sapply(function(x) x[2]) %>% sub(".smr", "", .),
    full_paths=list.files(dir, full.names = T, pattern = "smr$")
  )

  if (length(unique(df_tmp$exposures)) > 1) {
    cli::cli_alert_info("{format(Sys.time(), '%H:%M')}: Deal with {.emph {length(unique(exposures))}} exposure{?s} seperately.")
  }
  if (length(unique(df_tmp$outcomes)) > 1) {
    cli::cli_alert_info("Deal with {.emph {length(unique(outcomes))}} outcome{?s} seperately.")
  }
  if (out_dir == "") {
    out_dir <- file.path(dir, "chr_combined")
    cli::cli_alert_success("Out dir set to >>> {.path {out_dir}}")
  }
  if (!file.exists(out_dir)) {dir.create(out_dir)}

  for (exposure in unique(df_tmp$exposures)) {
    for (outcome in unique(df_tmp$outcomes)) {
      files_paths <- df_tmp %>% dplyr::filter(exposures == exposure, outcomes == outcome) %>% pull("full_paths")
      res <- lapply(files_paths, data.table::fread)
      res <- do.call(rbind, res)
      out_file <- paste0(out_dir, "/chr_combine_", exposure, "@", outcome, ".smr")
      cli::cli_alert_info("Write to file: {.path {out_file}}")
      data.table::fwrite(res, out_file)
    }
  }
  cli::cli_alert_success("All done!")
  return(NULL)
}

