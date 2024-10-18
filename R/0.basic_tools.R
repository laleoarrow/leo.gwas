#' Get ID
#'
#' Get ID using CHR, BP, A2 (REF/Non-effect), A1 (ALT/Effect)
#'
#' @keywords basic
#'
#' @param x A data.frame that must contain the columns CHR, BP, A2, and A1.
#' Each column represents:
#' - CHR: Chromosome number
#' - BP/POS: Base pair position
#' - A2: Reference allele (non-effect allele)
#' - A1: Alternative allele (effect allele)
#'
#' @return A df data.frame with ID column and the character counts of A1 and A2
#' @export
get_id <- function(x) {
  position_col <- if ("BP" %in% colnames(x)) {
    "BP"
  } else if ("POS" %in% colnames(x)) {
    "POS"
  } else {
    stop("Neither 'BP' nor 'POS' column found in the dataframe.")
  }
  x <- x %>%
    mutate(ID = paste(CHR, !!sym(position_col), A2, A1, sep = ":"),
           A1_n = nchar(A1),
           A2_n = nchar(A2))
  return(x)
}

#' Convert CHR:BP to rsID
#'
#' This function takes a dataframe that must include columns 'CHR' and 'BP',
#' and it appends the corresponding rsID by querying the SNPlocs.Hsapiens.dbSNP155.GRCh37 database.
#' The function returns a dataframe with the rsIDs included.
#'
#' @param df A dataframe containing at least the columns 'CHR' and 'BP'.
#' @param ref str indicating ref version.
#' - "GRCh37" for GRCh37 Ref panel
#' - "GRCh38" for GRCh38 Ref panel
#' @return A dataframe with an additional 'RefSNP_id' column which contains the rsIDs.
#' @importFrom data.table ":="
#' @examples
#' pacman::p_load(data.table, BSgenome, leo.gwas)
#' library("SNPlocs.Hsapiens.dbSNP155.GRCh37") # for GRCh37
#' library("SNPlocs.Hsapiens.dbSNP155.GRCh38") # for GRCh38
#' df <- data.frame(
#'   CHR = c(1, 1),
#'   BP = c(15211, 15820)
#' )
#' result <- add_rsid(df); result
#' @export
add_rsid <- function(dat, ref = "GRCh37") {
  if (!("CHR" %in% colnames(dat)) || !("BP" %in% colnames(dat))) {
    stop("DataFrame must contain 'CHR' and 'BP' columns")
  }

  # Convert dat to data.table if it is not one already
  # library(data.table)
  if (!data.table::is.data.table(dat)) {data.table::setDT(dat)}
  dat[, ranges := paste0(CHR, ":", BP, "-", BP)]

  # Load SNP data - assuming GRCh37, modify if using GRCh38
  if (ref == "GRCh37") {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37
  } else {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
  }

  # Find overlaps
  message(paste0("Translating RSID using: \n - BSgenome::snpsByOverlaps() \n - With ", snps@data_pkgname))
  snp.res <- BSgenome::snpsByOverlaps(snps, GRanges(dat$ranges))

  # Convert results to data.table
  snp.res.dt <- as.data.table(snp.res)
  snp.res.dt[, ranges := stringr::str_c(seqnames, ":", pos, "-", pos)]

  # Merge data
  trans.dat <- merge(snp.res.dt, dat, by = "ranges")
  columns_to_remove <- c("strand", "alleles_as_ambig", "ranges", "seqnames", "pos")
  trans.dat[, (columns_to_remove) := NULL]; gc()

  # Drop NA rsid and return
  # trans.dat <- trans.dat %>% drop_na(RefSNP_id)
  leo.gwas::leo_message("Remember to check if there is any NA in the RefSNP_id column.")
  leo.gwas::leo_message(">>> table(is.na(dat$RefSNP_id))")
  leo.gwas::leo_message(">>> vkh_meta %>% map_dbl(~sum(is.na(.)))")

  return(trans.dat)
}

#' Add Chromosome Position (CHR and POS) to Dataset Based on rsID
#'
#' This function takes a dataset with a column containing rsIDs (SNP IDs) and
#' adds the corresponding chromosome (CHR) and position (POS) information.
#' It queries the SNPlocs.Hsapiens.dbSNP155.GRCh37 database (or GRCh38 if specified) to retrieve the genomic positions.
#' The function returns a dataframe with the additional 'CHR' and 'POS' columns appended.
#'
#' @param dat A dataframe containing at least a column with SNP IDs (rsIDs).
#' @param snp_col A string indicating the column name containing SNP IDs (default is "SNP").
#' @param ref A string indicating the reference genome version. Default is "GRCh37", can also use "GRCh38".
#' @return A dataframe with additional 'CHR' and 'POS' columns.
#' @importFrom data.table ":="
#' @examples
#' pacman::p_load(data.table, dplyr, BSgenome, SNPlocs.Hsapiens.dbSNP155.GRCh37)
#' zuo_ref <- fread("/path/to/1KG-EAS-EAF.txt.gz") # Input dataset with rsID (SNP column)
#' result <- add_chrpos(zuo_ref, snp_col = "SNP", ref = "GRCh37")
#' @export
add_chrpos <- function(dat, snp_col = "SNP", ref = "GRCh37") {

  # Check if the SNP column exists in the dataframe
  if (!(snp_col %in% colnames(dat))) {
    stop(paste0("The column '", snp_col, "' is not found in the dataframe."))
  }

  # Load SNP reference data
  if (ref == "GRCh37") {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37
    leo.gwas::leo_message("ℹ Loading SNPlocs.Hsapiens.dbSNP155.GRCh37 database.")
  } else if (ref == "GRCh38") {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
    leo.gwas::leo_message("ℹ Loading SNPlocs.Hsapiens.dbSNP155.GRCh38 database.")
  } else {
    stop("Invalid reference version. Please choose 'GRCh37' or 'GRCh38'.")
  }

  # Get CHR and POS for the SNPs from the SNPlocs database
  snp_info <- data.table::setDT(data.frame(snpsById(snps, dat[[snp_col]], ifnotfound = "drop"))) %>%
    dplyr::select(seqnames, pos, RefSNP_id) %>%
    dplyr::rename(CHR = seqnames, POS = pos, SNP = RefSNP_id)

  # Merge the SNP info with the original dataset, keeping all original columns
  merged_data <- merge(dat, snp_info, by.x = snp_col, by.y = "SNP", all.x = TRUE) %>%
    select(SNP, CHR, POS, everything())

  # Check for missing SNPs (those not found in the SNP reference)
  missing_snps <- merged_data[is.na(merged_data$CHR), ]$SNP
  if (length(missing_snps) > 0) {
    message(paste0("Cannot find the following SNPs:\n", paste(missing_snps, collapse = ", "), "\n"))
  }

  # Return the final dataset with CHR and POS columns added
  return(merged_data)
}

#' Give precise p-value for chi-square test
#'
#' Sometimes you get a p-value of 0 when you perform a chi-square test or other analysis.
#' This is because the p-value is so small that R rounds it to 0. This function gives you
#' a more precise p-value.
#'
#' @param chisq_value numeric, chi-square value
#' @param df numeric, degree of freedom
#' @param digits numeric, digits for output to illustrate
#' @param prec numeric, precision for mpfr() function
#'
#' @return A precise p-value in scientific format
#' @export
#'
#' @examples
#' # install.packages("Rmpfr")
#' library(Rmpfr)
#' chisq_value <- 2629; df <- 1
#' p_value <- chisq_p_value(chisq_value, df)
#' print(p_value)
chisq_p_value <- function(chisq_value, df, digits = 4, prec = 100) {
  log_p_value <- pchisq(chisq_value, df, lower.tail = FALSE, log.p = TRUE)
  p_value <- exp(mpfr(log_p_value, prec = prec))
  # p_value <- format(p_value, scientific = TRUE, digits = digits)
  return(p_value)
}

#' Exclude HLA region from genomic summary data
#'
#' This function filters out entries within the HLA region on a specified chromosome
#' and position range. The default HLA region is set to chromosome 6, between 25Mb and 34Mb.
#' Custom chromosome and position bounds can be specified.
#'
#' @param data A data frame containing genomic data.
#' @param chromosome_col The name of the column representing chromosome numbers (default is CHR).
#' @param position_col The name of the column representing genomic positions (default is BP).
#' @param lower_bound The lower boundary of the HLA region in base pairs (default is 25e6).
#' @param upper_bound The upper boundary of the HLA region in base pairs (default is 34e6).
#'
#' @return A data frame excluding rows that fall within the specified HLA region.
#' @importFrom dplyr filter
#' @export
#' @examples
#' example_data <- data.frame(
#'   SNP = c("rs1", "rs2", "rs3", "rs4"),
#'   CHR = c(6, 6, 6, 7),
#'   POS = c(26000000, 33000000, 35000000, 29000000)
#' )
#' result_data <- exclude_HLA(example_data, chromosome_col="CHR", position_col="POS")
#' print(result_data)
exclude_HLA <- function(data, chromosome_col="CHR", position_col="BP", lower_bound=25e6, upper_bound=34e6) {
  filtered_data <- data %>%
    dplyr::filter(!(!!sym(chromosome_col) == 6 &
                      !!sym(position_col) >= lower_bound &
                      !!sym(position_col) <= upper_bound))
  return(filtered_data)
}
