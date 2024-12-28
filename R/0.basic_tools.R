#' Give Messages with my color
#'
#' @param ... The messages you wanna messgae, which will be pasted together
#' @param color prefered color
#' @param return Logical. If TRUE, returns the formatted string. If FALSE, prints directly.
#'
#' @export
#' @examples
#' leo_message("This is a pink message")
#' leo_message("This is a green message", "\nhaha", color = "32")
#' leo_message("This is a blue ", "message", color = "34")
#' leo_message("This is a light purple message", color = "95")
leo_message <- function(..., color = "31", return = FALSE) {
  message_content <- paste0(..., collapse = " ")
  formatted_message <- paste0("\033[", color, "m", message_content, "\033[0m")
  if (return) {
    return(formatted_message)
  } else {
    message(formatted_message)
  }
}

#' Log Messages with Timestamps
#'
#' Logs messages with timestamps
#' The messages are styled using the `cli` package for enhanced readability.
#' This function can not deal with {} function in the `cli` package.
#'
#' @param ... The message string to log, which will be pasted together.
#' @param level The log level. Options are `"info"`, `"success"`, `"warning"`, and `"danger"`.
#' @param levels All levels that is now supported.
#'
#' @return No return value. Outputs a formatted log message with a timestamp.
#' @export
#'
#' @examples
#' leo_log("Processing the", n1, "and", n2, "files.")
#' leo_log("Task completed successfully!", level = "success")
#' leo_log("Potential issue detected.", level = "warning")
#' leo_log("Error occurred during processing!", level = "danger")
leo_log <- function(..., level = "info", levels = c("info", "success", "warning", "danger")) {
  msg <- paste(..., collapse = " ");level <- match.arg(levels);timestamp <- format(Sys.time(), '%H:%M')
  colored_timestamp <- switch(level,
    info    = leo_message(paste0("[",timestamp,"]"), color = "36", return = TRUE), # blue
    success = leo_message(paste0("[",timestamp,"]"), color = "32", return = TRUE), # green
    warning = leo_message(paste0("[",timestamp,"]"), color = "33", return = TRUE), # red
    danger  = leo_message(paste0("[",timestamp,"]"), color = "31", return = TRUE)  # red
  )
  switch(level,
         info    = cli::cli_alert_info("{colored_timestamp} {msg}"),
         success = cli::cli_alert_success("{colored_timestamp} {msg}"),
         warning = cli::cli_alert_warning("{colored_timestamp} {msg}"),
         danger  = cli::cli_alert_danger("{colored_timestamp} {msg}")
  )
}


#' Get Unique Identifier for Genetic Data
#'
#' Get ID using CHR, BP, A2 (REF/Non-effect), A1 (ALT/Effect)
#'
#' @param x A data.frame that must contain the columns CHR, BP, A2, and A1.
#' Each column represents:
#' - CHR: Chromosome (It can be any in c("chrom", "CHR", "Chromosome", "chromosome", "Chr"))
#' - BP/POS: Base pair position (It can be any in c("pos", "POS", "position", "BP", "Position", "Bp"))
#' - A2: Reference allele/non-effect allele (It can be any in c("A2", "Allele2", "allele2", "a2", "REF", "Ref", "ref", "Non-effect"))
#' - A1: Alternative allele/effect allele (It can be any in c("A1", "Allele1", "allele1", "a1", "ALT", "Alt", "alt", "Effect"))
#' @param count_A1_A2 If T, will count the number of characters in A1 and A2
#'
#' @return A data.frame with an additional 'ID' column (if count_A1_A2=T, containing unique identifiers and character counts of A1 and A2)
#' @import dplyr
#' @examples
#' df <- data.frame(chrom = c(1, 1, 2), pos = c(12345, 54321, 11111), A1 = c("A", "T", "G"), A2 = c("G", "C", "A"))
#' get_id(df); get_id(df, count_A1_A2 = T)
#' @export
get_id <- function(x, count_A1_A2 = F) {
  # require(dplyr)
  # possible colnames for CHR and POS
  chrom_cols <- c("chrom", "CHR", "Chromosome", "chromosome", "Chr")
  pos_cols <- c("pos", "POS", "position", "BP", "Position", "Bp")
  a1_cols <- c("A1", "Allele1", "allele1", "a1", "ALT", "Alt", "alt", "Effect")
  a2_cols <- c("A2", "Allele2", "allele2", "a2", "REF", "Ref", "ref", "Non-effect")

  chrom_col <- intersect(chrom_cols, names(x)); if(length(chrom_col) == 0){stop("No chromosome column found in the dataframe.")}
  chrom_col <- chrom_col[1]

  pos_col <- intersect(pos_cols, names(x)); if(length(pos_col) == 0){stop("No position column found in the dataframe.")}
  pos_col <- pos_col[1]

  a1_col <- intersect(a1_cols, names(x)); if(length(a1_col) == 0){stop("No A1 column found in the dataframe.")}
  a1_col <- a1_col[1]

  a2_col <- intersect(a2_cols, names(x)); if(length(a2_col) == 0){stop("No A2 column found in the dataframe.")}
  a2_col <- a2_col[1]

  if (count_A1_A2 == TRUE) {
    x <- x %>% dplyr::mutate(ID = paste(.data[[chrom_col]], .data[[pos_col]], .data[[a2_col]], .data[[a1_col]], sep = ":"),
                             A1_n = base::nchar(.data[[a1_col]]), A2_n = nchar(.data[[a2_col]]))
  } else {
    x <- x %>% dplyr::mutate(ID = paste(.data[[chrom_col]], .data[[pos_col]], .data[[a2_col]], .data[[a1_col]], sep = ":"))
  }

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

#' Count or Identify Matches of a Pattern in a Vector
#'
#' This function counts the number of elements in a vector that contain a given pattern,
#' or returns a logical vector indicating which elements match the pattern.
#'
#' @param vector A character vector to be searched.
#' @param pattern The pattern to match (regular expression).
#' @param return A string specifying the return type: "count" for number of matches, or "logi" for a logical vector.
#'
#' @return An integer representing the number of matches if return is "count",
#'         or a logical vector indicating which elements match the pattern if return is "logi".
#' @export
#' @examples
#' vec <- c("abc", "def", "xyz", "abcd")
#' count_matching_elements(vec, "abc", return = "count")  # Returns 2
#' count_matching_elements(vec, "abc", return = "logi")   # Returns c(TRUE, FALSE, FALSE, TRUE)
count_matching_elements <- function(vector, pattern, return = "count") {
  # Ensure the input is a character vector
  if (!is.character(vector)) {
    stop("Input must be a character vector.")
  }

  # Create a logical vector for matching elements
  matches <- grepl(pattern, vector, ignore.case = TRUE)

  if (return == "count") {
    return(sum(matches))  # Return the count of matching elements
  } else if (return == "logi") {
    return(matches)  # Return the logical vector
  } else {
    stop("Invalid return value. Use 'count' or 'logi'.")
  }
}

#' Count or Identify Duplicates in a Vector
#'
#' This function counts the number of duplicate elements in a vector,
#' or returns a logical vector indicating which elements are duplicates.
#'
#' @param vector A vector to be checked for duplicates.
#' @param return A string specifying the return type: "count" for number of duplicates, or "logi" for a logical vector.
#'
#' @return An integer representing the number of duplicates if return is "count",
#'         or a logical vector indicating which elements are duplicates if return is "logi".
#' @export
#' @examples
#' vec <- c("a", "b", "c", "a", "b", "d")
#' count_duplicate_element(vec, return = "count")  # Returns 4
#' count_duplicate_element(vec, return = "logi")   # Returns c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
count_duplicate_element <- function(vector, return = "count") {
  # Find duplicates
  duplicate_indices <- duplicated(vector) | duplicated(vector, fromLast = TRUE)

  if (return == "count") {
    return(sum(duplicate_indices))  # Return the count of duplicates
  } else if (return == "logi") {
    return(duplicate_indices)  # Return the logical vector
  } else {
    stop("Invalid return value. Use 'count' or 'logi'.")
  }
}


#' Across a df to count na
#'
#' This function summarize NA values in each column of a data frame.
#'
#' @param df a data frame
#' @return a data frame with the number of NA values in each column
#' @export
#' @examples
#' df <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 2, 3, 4))
#' summarize_na(df)
across_df_na <- function(df){
  df %>% dplyr::summarise(across(everything(), ~sum(is.na(.))))
}

#' Across a df to count TRUE and FALSE
#'
#' This function summarizes both TRUE and FALSE values in each column of a data frame.
#'
#' @param df a data frame
#' @param type "T" for TRUE counts (default), "F" for FALSE counts
#' @return a data frame with the number of TRUE or FALSE values in each column
#' @export
#' @examples
#' df <- data.frame(a = c(TRUE, FALSE, TRUE, TRUE), b = c(FALSE, TRUE, TRUE, TRUE))
#' across_df_TF(df) # Count TRUE (default)
#' across_df_TF(df, "F") # Count FALSE
across_df_TF <- function(df, type = "T"){
  type <- match.arg(type, c("T", "F"))
  if(type == "T"){
    df %>% dplyr::summarise(across(everything(), ~sum(.)))
  } else {
    df %>% dplyr::summarise(across(everything(), ~sum(!.)))
  }
}
