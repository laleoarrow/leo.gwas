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
