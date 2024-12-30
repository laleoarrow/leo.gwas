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


# ---- Below is how you aggregate all SMR results for XWAS ----
# smr result should be all stored in a dir in the following format:
# - level 1: qtl type (e.g. mqtl, eqtl, ...)
# - level 2: qtl specific source (e.g. GTEx49, ...)

# 0. Combine SMR res of all CHR from one source ----
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
  leo_log("Combine SMR results for all chromosomes.")
  df_tmp <- dplyr::tibble(
    files=list.files(dir, full.names = F, pattern = "smr$"),
    exposures=strsplit(files, "@") %>% sapply(function(x) x[1]) %>% sub("_chr[0-9]+", "", .),
    outcomes=strsplit(files, "@") %>% sapply(function(x) x[2]) %>% sub(".smr", "", .),
    full_paths=list.files(dir, full.names = T, pattern = "smr$")
  )

  if (length(unique(df_tmp$exposures)) > 1) {
    cli::cli_alert_info("Deal with {.emph {length(unique(df_tmp$exposures))}} exposure{?s} seperately.")
  }
  if (length(unique(df_tmp$outcomes)) > 1) {
    cli::cli_alert_info("Deal with {.emph {length(unique(df_tmp$outcomes))}} outcome{?s} seperately.")
  }
  if (out_dir == "") {
    out_dir <- file.path(dir, "chr_combined")
    leo_log("Out dir set to >>>", out_dir, level = "success")
  }
  if (!file.exists(out_dir)) {dir.create(out_dir)}

  for (exposure in unique(df_tmp$exposures)) {
    for (outcome in unique(df_tmp$outcomes)) {
      files_paths <- df_tmp %>% dplyr::filter(exposures == exposure, outcomes == outcome) %>% pull("full_paths")
      res <- lapply(files_paths, data.table::fread)
      res <- do.call(rbind, res)
      out_file <- paste0(out_dir, "/chr_combine_", exposure, "@", outcome, ".smr")
      cli::cli_alert_info("Write to file: {.path {out_file}}")
      data.table::fwrite(res, out_file, sep = "\t")
    }
  }
  leo_log("ALL DONE!", level = "success")
  return(NULL)
}

# 1. FDR/Bonferroni ----
#' Adjust SMR Results with FDR and Bonferroni Corrections
#'
#' It applies FDR/Bonferroni corrections for a single SMR result.
#' The corrected results are saved in an `fdr` subdirectory within the output directory.
#'
#' @param smr_result_path Character. The path containing SMR result file. Files should have the `.smr` extension.
#' @param out_dir Character. The output directory where adjusted files will be saved.
#'                If not specified, defaults to creat a dir named `fdr` within dirname(smr_result_path).
#' @param writePath
#' @param write_out Logical. Whether to write the adjusted results to a file. Default is `TRUE`.
#' @return NULL
#' @examples
#' \dontrun{
#' leo_smr_adjust("~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined/chr_combine_sQTL_Adipose_Subcutaneous@iri3.smr",
#'writePath = "",
#'out_dir = "~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49")
#'leo_smr_adjust("~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined/chr_combine_sQTL_Adipose_Subcutaneous@iri3.smr",
#'               writePath = "./haha.fdr",
#'               out_dir = "")
#' }
#' @importFrom cli cli_alert_info
#' @importFrom vroom vroom vroom_write
#' @importFrom dplyr mutate %>%
#' @importFrom stringr str_replace
#' @export
leo_smr_adjust <- function(smr_result_path, writePath = "", out_dir = "", write_out = T) {
  smr_result <- vroom::vroom(smr_result_path, show_col_types = F) %>%
    mutate(HLA_Probe = ifelse( (ProbeChr == 6 & Probe_bp >= 25000000 & Probe_bp <= 34000000), "Yes", "No") )
  nrow_HLA_probe <- smr_result %>% filter(HLA_Probe == "Yes") %>% nrow()
  cli::cli_alert_info("Filtering out {.emph {nrow_HLA_probe}} probe{?s} for {basename(smr_result_path)}")
  smr_result <- smr_result %>% filter(HLA_Probe == "No") %>%
    mutate(Pass_HEIDI = ifelse(p_HEIDI>=0.05, "Pass", "Fail"),
           N_probe = nrow(.),
           FDR = p.adjust(p_SMR, method = "BH"),
           Pass_FDR = ifelse(FDR < 0.05, "Pass", "Fail"),
           Bonferroni = p.adjust(p_SMR, method = "bonferroni"),
           Pass_Bonferroni = ifelse(Bonferroni < 0.05, "Pass", "Fail"))

  if (write_out) {
    basename <- basename(smr_result_path) %>% gsub(".smr", ".fdr", .)
    dirname <- dirname(smr_result_path)

    if (writePath == "") {
      if (out_dir == "") {
        cli::cli_alert_warning("No writePath & out_dir set")
        out_dir <- file.path(dirname, "fdr")
      }
      if (!file.exists(out_dir)) {dir.create(out_dir)}
      writePath <- file.path(out_dir, basename)
    } else {
      writePath <- writePath # this will overide the basename setting
    }

    cli::cli_alert_success(" - Writing to >>> {.path {writePath}}")
    vroom::vroom_write(smr_result, writePath, delim = "\t")
    return(NULL)
  } else {
    return(smr_result)
  }
}


#' Batch Adjust SMR Results with FDR and Bonferroni Corrections
#'
#' This function applies FDR/Bonferroni corrections to all SMR results within a specified directory.
#'
#' @param dir Character. The directory containing SMR result files. Files should have the `.smr` extension and follow the naming convention `exposure@outcome.smr`.
#' @param out_dir Character. The output directory where the adjusted SMR files will be saved.
#'                If not specified, defaults to creating a `fdr` directory within the input directory.
#' @param ... Additional arguments to be passed to `leo_smr_adjust`.
#'
#' @return NULL.
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#' @importFrom dplyr filter pull
#' @export
leo_smr_adjust_loop <- function(dir, out_dir="", ...) {
  leo_log("Adjusting SMR results for all files in the >>>", dir, "<<< directory.")
  if (out_dir == "") {
    dirname <- dirname(dir)
    out_dir <- file.path(dirname, "fdr")
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    cli::cli_alert_warning("out_dir set to {out_dir} as there are no pre-set out_dir.")
  } else {
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir = out_dir
  }

  # list all smr file
  smr_files <- list.files(path = dir, pattern = "\\.smr$", full.names = TRUE) # Get list of .smr files in the folder
  if (length(smr_files) == 0) {
    leo_log("No `.smr` files found in the directory.", level = "danger")
    return(NULL)
  }
  for (smr_file in smr_files) { # Loop through each .smr file and apply leo_smr_adjust function
    cli::cli_alert_info(" - Processing: {smr_file}")
    leo_smr_adjust(smr_file, out_dir = out_dir, ...)
  }
  leo_log("ALL DONE!", level = "success")
}
# leo_smr_adjust_loop(dir     = "~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined",
#                     out_dir = "~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/fdr")
