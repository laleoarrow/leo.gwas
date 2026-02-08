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
#' @importFrom dplyr mutate filter select pull group_by summarise left_join inner_join everything .data
#' @importFrom stats complete.cases cor.test median p.adjust pchisq t.test
#' @importFrom grDevices dev.off pdf
#' @examples
#' df <- data.frame(chrom = c(1, 1, 2), pos = c(12345, 54321, 11111),
#'                  A1 = c("A", "T", "G"), A2 = c("G", "C", "A"))
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

#' Leo batch iterator builder
#'
#' Sometimes you need to handle a large number of iterations, but multi-core
#' parallel computing can be tricky—memory usage may grow beyond what is
#' actually required. This function helps by batching the iterations, so you
#' only process a limited number of elements at a time, preventing excessive
#' RAM consumption.
#'
#' @param elements A vector or list to be iterated.
#' @param batch_size Number of elements per batch.
#' @details Each call to the returned function yields the next batch. Returns \code{NULL} when no more elements remain.
#' @return A function (no arguments). Repeated calls produce successive batches or \code{NULL} if finished.
#' @examples
#' # Suppose you have 25 elements and want to batch them in groups of 6
#' nums <- 1:25
#' it <- leo_iterator(nums, 6)
#' while (TRUE) {
#'   batch <- it()
#'   if (is.null(batch)) break
#'
#'   # add your parallel process steps
#'   print(batch)
#' }
#' @export
leo_iterator <- function(elements, batch_size) {
  total_elements <- length(elements)
  total_rounds <- ceiling(total_elements / batch_size)
  current_round <- 0

  leo_log(
    "Creating batch iterator with", total_elements, "elements in total; ",
    "batch size =", batch_size, "; total rounds =", total_rounds, "\n"
  )

  function() {
    # If we've already returned all possible batches, return NULL
    if (current_round >= total_rounds) {
      return(NULL)
    }
    current_round <<- current_round + 1
    start_index <- (current_round - 1) * batch_size + 1
    end_index   <- min(total_elements, start_index + batch_size - 1)

    # Print iteration info
    leo_log("Round", current_round, "of", total_rounds,
            "| elements from index", start_index, "to", end_index, "\n",
            level = "success")

    elements[start_index:end_index]
  }
}

# ---------------------- basic calculation ----------------------
#' Calculate Correlation between Two Vectors
#'
#' This function calculates the Spearman (default) or Pearson correlation coefficient
#' and its associated p-value between two vectors.
#' It automatically handles missing values.
#'
#' @param vector_x A numeric vector.
#' @param vector_y A numeric vector of the same length as \code{vector_x}.
#' @param method A character string specifying the correlation method ("spearman" or "pearson").
#'               Defaults to "spearman".
#' @param ... Pass to \code{\link[stats]{cor.test}}.
#'
#' @return A list with the correlation coefficient and p-value.
#' @export
#' @seealso \code{\link{correlation_draw}} for plotting correlation results.
#' @examples
#' vector_x <- c(1, 2, 2, 4, 5)
#' vector_y <- c(5, 6, 7, 8, 7)
#' result <- correlation_calculate(vector_x, vector_y, method = "pearson")
correlation_calculate <- function(vector_x, vector_y, method = "spearman", ...) {
  if (!method %in% c("spearman", "pearson", "kendall")) {
    stop("Invalid method. Choose 'spearman', 'pearson', or 'kendall'.")
  }

  # Remove NA values from both vectors consistently
  valid_indices <- complete.cases(vector_x, vector_y)
  vector_x <- vector_x[valid_indices]
  vector_y <- vector_y[valid_indices]

  # Use cor.test to get correlation coefficient and p-value
  correlation_test <- cor.test(vector_x, vector_y, method = method, use = "complete.obs", ...)

  # Extract correlation coefficient and p-value
  result <- data.frame(
    correlation_coefficient = as.numeric(correlation_test$estimate),
    p_value = correlation_test$p.value
  )

  return(result)
}
