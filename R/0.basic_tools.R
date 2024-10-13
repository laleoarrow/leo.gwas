#' Get ID
#'
#' Get ID using CHR, BP, A2 (REF/Non-effect), A1 (ALT/Effect)
#'
#' @keywords basic
#'
#' @param x A data.frame that must contain the columns CHR, BP, A2, and A1.
#' Each column represents:
#' - CHR: Chromosome number
#' - BP: Base pair position
#' - A2: Reference allele (non-effect allele)
#' - A1: Alternative allele (effect allele)
#'
#' @return A df data.frame with ID column and the character counts of A1 and A2
#' @export
get_id <- function(x) {
  x <- x %>%
    mutate(ID = paste(CHR, BP, A2, A1, sep = ":"),
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
