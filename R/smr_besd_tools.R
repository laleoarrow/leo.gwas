# This contains need-to-have functions for making BESD files for SMR and HEIDI



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
#'                If not specified, defaults to create a dir named `fdr` within dir-name (smr_result_path).
#' @param writePath Character. The path to write the adjusted results.
#' @param hla Logical. Whether to pre-exclude the HLA region probes. Default is `False`.
#' @param write_out Logical. Whether to write the adjusted results to a file. Default is `TRUE`.
#'
#' @return NULL
#' @examples
#' \dontrun{
#' leo_smr_adjust("~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined/chr_combine_sQTL_Adipose_Subcutaneous@iri3.smr",
#'                writePath = "", out_dir = "~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49")
#' leo_smr_adjust("~/project/iridocyclitis/output/smr-t2d/sqtl/GTEx49/chr_combined/chr_combine_sQTL_Adipose_Subcutaneous@iri3.smr",
#'                writePath = "./haha.fdr", out_dir = "")
#' }
#' @importFrom cli cli_alert_info
#' @importFrom vroom vroom vroom_write
#' @importFrom dplyr mutate %>% filter
#' @importFrom stringr str_replace
#' @export
leo_smr_adjust <- function(smr_result_path, writePath = "", out_dir = "", hla = F, write_out = T) {
  smr_result <- vroom::vroom(smr_result_path, show_col_types = F)
  if (hla) {
    smr_result <- smr_result %>%
      dplyr::mutate(HLA_Probe = ifelse((ProbeChr == 6 & Probe_bp >= 25000000 & Probe_bp <= 34000000), "Yes", "No")) %>%
      dplyr::filter(HLA_Probe == "No")
    nrow_HLA_probe <- smr_result %>% filter(HLA_Probe == "Yes") %>% nrow()
    cli::cli_alert_info(" - Filtering out {.emph {nrow_HLA_probe}} probe{?s} for {basename(smr_result_path)}") # seems unnecessary
  }
  smr_result <- smr_result %>%
    mutate(Pass_HEIDI = ifelse(p_HEIDI >= 0.05, "Pass", "Fail"),
           N_probe = nrow(.),
           FDR = p.adjust(p_SMR, method = "BH", n = N_probe[1]),
           Pass_FDR = ifelse(FDR < 0.05, "Pass", "Fail"),
           Bonferroni = p.adjust(p_SMR, method = "bonferroni", n = N_probe[1]),
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
    return(invisible(NULL))
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
#' @examples
#' \dontrun{
#' leo_smr_adjust_loop(dir     = "~/project/iridocyclitis/output/smr/sqtl/GTEx49/chr_combined",
#'                     out_dir = "~/project/iridocyclitis/output/smr/sqtl/GTEx49/fdr")
#' }
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#' @importFrom dplyr filter pull
#' @export
leo_smr_adjust_loop <- function(dir, out_dir="", ...) {
  leo_log("Adjusting SMR results for all files in the >>>", dir, "<<< directory.")
  if (out_dir == "") {
    dirname <- dirname(dir)
    out_dir <- file.path(dirname, "fdr")
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    cli::cli_alert_warning("out_dir set to {.path {out_dir}} as there are no pre-set out_dir.")
  } else {
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir = out_dir
  }

  # list all smr file
  smr_files <- list.files(path = dir, pattern = "\\.smr$", full.names = TRUE) # Get list of .smr files in the folder
  if (length(smr_files) == 0) {
    leo_log("No `.smr` files found in the directory.", level = "danger")
    return(invisible(NULL))
  }
  for (smr_file in smr_files) { # Loop through each .smr file and apply leo_smr_adjust function
    cli::cli_alert_info(" - Processing: {smr_file}")
    leo_smr_adjust(smr_file, out_dir = out_dir, ...)
  }
  leo_log("ALL DONE!", level = "success")
  return(invisible(NULL))
}

# 2. Combine all results for each outcomes ----
#' Combine `.fdr` files for one or multiple outcomes (seperately)
#'
#' @param dir       Character. The parent folder that contains subfolders with fdr files.
#' @param outcome   Character or character vector. One or more outcomes to search in file names.
#' @param out_dir   Character. Output folder; default is to create "combine_1outcome" under \code{dir}.
#' @return          NULL. This function writes combined files to disk directly.
#' @examples
#' \dontrun{
#' combine_smr_res_1outcome("/Users/leoarrow/project/iridocyclitis/output/smr", "iri3")
#' }
#' @importFrom cli cli_alert_danger cli_alert_success cli_alert_warning cli_alert_info
#' @importFrom dplyr bind_rows select everything
#' @importFrom data.table fwrite
#' @importFrom vroom vroom
#' @export
combine_smr_res_1outcome <- function(dir, outcome, out_dir = file.path(dir, "combine_1outcome")) {
  date_str <- format(Sys.Date(), "%Y%m%d")
  cli::cat_rule(paste0("[",date_str,"] Combining SMR results"), col = "blue")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE) # Create out_dir if not exist

  # Find all .fdr files
  fdr_files <- list.files(
    path         = dir,
    pattern      = "\\.fdr$",
    recursive    = TRUE,
    full.names   = TRUE
  )
  if (length(fdr_files) == 0) {
    leo_log("No `.fdr` files found in the given directory", level = "danger")
    return(invisible(NULL))
  }
  # check if given oc all in the file.names
  extract_outcome_from_filename <- function(fp) {
    fn_no_ext <- basename(fp) %>% sub("\\.fdr$", "", .) %>% sub("\\.ma$", "", .)
    if (grepl("@", fn_no_ext)) {
      return(sub(".*@", "", fn_no_ext))
    } else if (grepl("_", fn_no_ext)) {
      return(sub(".*_", "", fn_no_ext))
    } else {
      return(fn_no_ext)
    }
  }
  file.outcomes <- unique(sapply(fdr_files, extract_outcome_from_filename))
  missing_outcomes <- setdiff(outcome, file.outcomes)
  if (length(missing_outcomes) > 0) {
    msg <- paste0("The following outcome(s) are not found: ", paste(missing_outcomes, collapse = ", "), " Continue?")
    if (interactive()) {
      proceed <- utils::askYesNo(
        msg, default = T, prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel")))
      )
      if (!isTRUE(proceed)) {
        cli::cli_alert_danger("Operation aborted by the user.")
        return(invisible(NULL))
      }
    } else {
      cli::cli_alert_warning("Missing outcomes: {paste(missing_outcomes, collapse = ', ')}. Skipping them.")
    }
  }

  # Loop each outcome
  for (oc in outcome) {
    cli::cli_alert_info("Processing outcome: {oc}")

    # Select files with this outcome
    sub_files <- fdr_files[grepl(oc, basename(fdr_files))]
    if (length(sub_files) == 0) {
      leo_log(" - No `.fdr` file matched: ", oc, level = "warning")
      next
    }

    # Read and combine
    df_list <- lapply(sub_files, function(fp) {
      tmp <- vroom::vroom(fp, show_col_types = F, progress = F) %>% dplyr::mutate(
        source_type = basename(dirname(dirname(dirname(fp)))),
        source      = basename(dirname(dirname(fp))),
        filename    = basename(fp))
    }) %>% as.list()
    merged_df <- do.call(dplyr::bind_rows, df_list) %>% dplyr::select(filename, source, source_type, dplyr::everything())

    # Output file
    out_file <- file.path(out_dir, paste0("smr_res_", date_str, "_", oc, ".all"))
    cli::cli_alert_success("Writing => {out_file}")
    data.table::fwrite(merged_df, file = out_file, sep = "\t")
  }
  leo_log("All Done!", level = "success")
  return(invisible(NULL))
}

# 3. Extract significant results from `.all` files ----
#' Extract significant results from `.all` files
#'
#' This function scans the `combine_1outcome` folder for `.all` files,
#' filters rows where \code{Pass_FDR == "Pass"} or \code{Pass_Bonferroni == "Pass"}
#' (depending on \code{pass_type} argument), and writes the significant subset
#' to a new file (e.g., `smr_res_YYYYMMDD_iri3.sig.all`).
#'
#' @param dir        Character. The main `combine_1outcome` folder from \func{combine_smr_res_1outcome}.
#' @param pass_type  Character. Which significance criterion to use:
#'                   one of \code{c("FDR", "Bonferroni", "both")}.
#'                   - "FDR": keep rows where \code{Pass_FDR == "Pass"}
#'                   - "Bonferroni": keep rows where \code{Pass_Bonferroni == "Pass"}
#'                   - "both": keep rows where either FDR or Bonferroni is "Pass"
#' @param out_dir    Character. Where to write significant results.
#'                   Default \code{dir/combine_1outcome/sig}.
#'
#' @return NULL. Writes `.sig.all` files to disk.
#'
#' @examples
#' \dontrun{
#' # For the directory "/Users/leoarrow/project/iridocyclitis/output/smr-t2d",
#' # after running combine_smr_res_1outcome, we have .all files in
#' # "smr-t2d/combine_1outcome".
#'
#' leo_smr_extract_sig_res(
#'   dir       = "/Users/leoarrow/project/iridocyclitis/output/smr-t2d",
#'   pass_type = "FDR"
#' )
#' }
#' @importFrom cli cli_alert_success cli_alert_info cli_alert_warning
#' @importFrom vroom vroom vroom_write
#' @importFrom dplyr filter
#' @export
leo_smr_extract_sig_res <- function(dir, pass_type = c("FDR", "Bonferroni", "both"), out_dir   = "") {
  cli::cat_rule("Extract significant results from `.all` files", col = "blue")
  pass_type <- match.arg(pass_type)
  if (!dir.exists(dir)) {
    cli::cli_alert_warning("No 'combine_1outcome' folder found in {dir}. Abort.")
    return(invisible(NULL))
  }
  # Default out_dir => dir/combine_1outcome/sig
  if (out_dir == "") {
    out_dir <- file.path(dir, "sig")
    cli::cli_alert_info("No out_dir specified. So out_dir set to {.path {out_dir}}")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # List all .all files, e.g. "smr_res_20241230_iri3.all"
  all_files <- list.files(
    path       = dir,
    pattern    = "\\.all$",
    full.names = TRUE
  )
  if (length(all_files) == 0) {
    cli::cli_alert_warning("No '.all' files found in {dir}")
    return(invisible(NULL))
  }

  # Loop through each .all file
  for (fpath in all_files) {
    fn <- basename(fpath)
    cli::cli_alert_info("Processing file: {fn}")
    df <- vroom::vroom(fpath, show_col_types = FALSE, progress = FALSE)

    # Filter by pass_type
    # if (!all(c("Pass_FDR", "Pass_Bonferroni") %in% names(df))) {
    #   cli::cli_alert_warning("File {fn} lacks 'Pass_FDR' or 'Pass_Bonferroni' columns. Skipped.")
    #   next
    # }
    df_sig <- switch(
      pass_type,
      "FDR" = dplyr::filter(df, .data$Pass_FDR == "Pass"),
      "Bonferroni" = dplyr::filter(df, .data$Pass_Bonferroni == "Pass"),
      "both" = dplyr::filter(df, .data$Pass_FDR == "Pass" | .data$Pass_Bonferroni == "Pass")
    ) %>% dplyr::filter(Pass_HEIDI == "Pass")

    if (nrow(df_sig) == 0) {
      cli::cli_alert_warning("No significant rows found in {fn} under '{pass_type}'. Skipped.")
      next
    }

    # Build & write output filename: e.g. "smr_res_20241230_iri3.sig.all"
    out_fn <- sub("\\.all$", ".sig.all", fn)
    out_fpath <- file.path(out_dir, out_fn)
    vroom::vroom_write(df_sig, out_fpath, delim = "\t")
    cli::cli_alert_success("Wrote significant results => {out_fpath} with {nrow(df_sig)} rows.")
  }
  return(invisible(NULL))
}

# 4. Shared ones ----
# In this section, we locate the shared XWAS results for different outcomes
#' Merge multiple SMR files and keep only shared probes
#'
#' This function finds intersected \code{probeID} across all provided files
#' and merges them into a single data frame using \code{inner_join}.
#'
#' @param dir        Character. A directory containing SMR files (e.g. ".sig.all" files).
#'                   If provided, the function will read all matching files in this directory.
#' @param file_paths Character vector. A set of absolute paths to SMR files. (Priority over \code{dir})
#'                   If this is non-NULL, \code{dir} is ignored.
#' @param pattern    Character. Regex pattern for searching files in \code{dir}. Default is "\\.sig\\.all$".
#' @param out_file   Character. If not empty, write the merged result to this path. Otherwise only return in R.
#'
#' @return A data frame (tibble) with shared probes. Also writes to \code{out_file} if \code{out_file != ""}.
#'
#' @examples
#' \dontrun{
#' # 1) Merge all .sig.all files in a folder:
#' merged_df <- leo_smr_merge_shared_probes(
#'   dir = "/Users/leoarrow/project/iridocyclitis/output/smr-t2d/combine_1outcome/sig"
#' )
#'
#' # 2) Merge specific files by absolute paths:
#' files_vec <- c(
#'   "/Users/leoarrow/project/iridocyclitis/output/smr/combine_1outcome/smr_res_20241230_iri3.all",
#'   "/Users/leoarrow/project/iridocyclitis/output/smr/combine_1outcome/smr_res_20241230_t1d1.all"
#' )
#' merged_df <- leo_smr_merge_shared_probes(file_paths = files_vec, out_file = "shared_probes_merged.tsv")
#' }
#'
#' @importFrom vroom vroom_write vroom
#' @importFrom dplyr inner_join group_by summarise filter select arrange
#' @importFrom purrr reduce
#' @export
leo_smr_merge_shared_probes <- function(
    dir         = NULL, file_paths  = NULL,
    pattern     = "\\.sig\\.all$",
    out_file    = "") {
  cli::cat_rule("Merge multiple SMR files and keep only shared probes", col = "blue")

  # 1) Identify files
  if (!is.null(file_paths)) {
    # If file_paths is provided, ignore dir
    leo_log("Reading files from 'file_paths': ", file_paths)
    smr_files <- file_paths
  } else {
    # Otherwise, search dir for pattern
    if (is.null(dir) || !dir.exists(dir)) stop("Please provide a valid 'dir' or a non-empty 'file_paths'.")
    leo_log("Reading files in 'dir': ", dir)
    smr_files <- list.files(
      path       = dir,
      pattern    = pattern,
      full.names = TRUE
    )
  }
  if (length(smr_files) < 2) stop("At least 2 files are required to find shared probes.")

  # 2) Read all files
  df_list <- lapply(smr_files, function(fp) {
    vroom::vroom(fp, show_col_types = FALSE, progress = FALSE)
  })

  # 3) Iteratively merge with inner_join by "probeID"
  #    This keeps only the probeID rows that appear in *all* files.
  #    If some file doesn't have "probeID" column, it will fail.
  merged_df <- purrr::reduce(
    .x  = df_list,
    .f  = function(x, y) dplyr::inner_join(x, y, by = c("source", "source_type", "probeID"))
  )

  # 4) Optionally write to out_file
  if (out_file != "") {
    vroom::vroom_write(merged_df, out_file, delim = "\t")
    cli::cli_alert_success("Merged data written to {out_file} with {nrow(merged_df)} rows.")
    return(invisible(NULL))
  } else {
    cli::cli_alert_info("Merged data has {nrow(merged_df)} shared rows across all files.")
    return(merged_df)
  }
}
# merged_df <- leo_smr_merge_shared_probes(
#     dir = "/Users/leoarrow/project/iridocyclitis/output/smr/combine_1outcome/sig"
#   )
# leo_smr_merge_shared_probes(
#   file_paths = c(
#     "/Users/leoarrow/project/iridocyclitis/output/smr/combine_1outcome/smr_res_20241230_iri3.all",
#     "/Users/leoarrow/project/iridocyclitis/output/smr-t2d/combine_1outcome/smr_res_20241230_iri3.all"
#   )
# )
