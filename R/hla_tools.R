# conditional analysis ----
# Global variables for R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(
  "SNP", "P", "CHISQ"
))

#' locate the significant SNP for conditional analysis
#'
#' @param x data.frame of the SNP information
#' @param environment.list tibble of the environment list (see the example for usage)
#' @param significance_level numeric, significance level, default is 5e-8
#'
#' @return message that informs the user of the independant SNP for the following conditional analysis
#' @export
#' @examples
#' \dontrun{
#' con_dir <- "/Users/leoarrow/project/VKH2024/data/zuo/con_su" 
#' # specify the directory to store the HLA original data and subsequent conditional analysis results data
#' files <- list.files(con_dir,full.names = T) \%>\% as.vector(); files # update it each time
#' x1 <- fread(files[1]) \%>\% arrange(desc(CHISQ)); head(x1) # ! for the first one, just read the data 
#' # and sort it by CHISQ/P value
#' x2 <- fread(files[2]) \%>\% arrange(P) # repeat it until no more independent signal can be found.
#' x3 <- fread(files[3]) \%>\% arrange(P)
#' x4 <- fread(files[4]) \%>\% arrange(P)
#' x5 <- fread(files[5]) \%>\% arrange(P)
#' x6 <- fread(files[6]) \%>\% arrange(P)
#' x7 <- fread(files[7]) \%>\% arrange(P)
#' x8 <- fread(files[8]) \%>\% arrange(P)
#'
#' env <- ls() # get the environment; this line if were put in main func will lead to error.
#' # load the environment
#' environment.list <- tibble(item = as.vector(grep("^x[0-9]+$", x = env, value = TRUE))) 
#' check_significant_SNP(x8, environment.list, significance_level = 5e-8)
#'
#' # You can mannually check the p-value of one SNP in previous environment.list
#' }
#' 
#'

check_significant_SNP <- function(x, environment.list, significance_level = 5e-8) {
  if (nrow(environment.list) == 0) return(message("No matching variables found in the environment."))

  for (i in 1:nrow(x)) { # i = 1
    con_SNP <- x$SNP[i]
    con_SNP_p <- x$P[i]

    if (con_SNP_p > significance_level) {
      message(paste("SNP: ", con_SNP, " is not significant already."))
      leo_message("The conditional analysis can be over now. So long!")
      break
    }

    message(paste("Locating in previous environment for <<<<<<<", con_SNP, ">>>>>>>"))
    all_pass <- lapply(1:(nrow(environment.list)-1), function(j){
      tmp_ref <- get(environment.list$item[j])
      tmp_ref_p <- tmp_ref %>% filter(SNP == con_SNP) %>% pull(P) %>% as.numeric()
      return(tmp_ref_p < significance_level)
    }) %>% do.call(all, .)

    if (all_pass) {
      message(paste0(" <<<< Passed for SNP: ", con_SNP, " >>>>"))
      leo_message("It can be used for the next SNP for conditional analysis")
      leo_log("A total of {i} search has been made.", level = "success")
      break
    } else {
      leo_message("Not this one! Let's move on to the next one.")
    }

    if (i == nrow(x)) {
      leo_message("The conditional analysis can be over now. So long!")
    }
  }
}

# others ----
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
#' @importFrom dplyr filter %>%
#' @export
#' @examples
#' example_data <- data.frame(
#'   SNP = c("rs1", "rs2", "rs3", "rs4"),
#'   CHR = c(6, 6, 6, 7),
#'   POS = c(26000000, 33000000, 35000000, 29000000)
#' )
#' result_data <- exclude_HLA(example_data, chromosome_col="CHR", position_col="POS")
#' print(result_data)
HLA_exclude <- function(data, chromosome_col="CHR", position_col="BP", lower_bound=25e6, upper_bound=34e6) {
  filtered_data <- data %>%
    dplyr::filter(!(!!sym(chromosome_col) == 6 &
                      !!sym(position_col) >= lower_bound &
                      !!sym(position_col) <= upper_bound))
  return(filtered_data)
}
#' Extract HLA region from genomic summary data
#'
#' This function extracts entries within the HLA region on a specified chromosome
#' and position range. The default HLA region is chromosome 6, between 25Mb and 34Mb.
#' Users may provide custom chromosome and region boundaries.
#'
#' @param data A data frame containing genomic data.
#' @param chromosome_col The name of the column representing chromosome numbers (default is CHR).
#' @param position_col The name of the column representing genomic positions (default is BP).
#' @param lower_bound The lower boundary of the HLA region in base pairs (default is 25e6).
#' @param upper_bound The upper boundary of the HLA region in base pairs (default is 34e6).
#'
#' @return A data frame including only rows within the specified HLA region.
#' @importFrom dplyr filter %>%
#' @export
#' @examples
#' example_data <- data.frame(
#'   SNP = c("rs1", "rs2", "rs3", "rs4"),
#'   CHR = c(6, 6, 6, 7),
#'   POS = c(26000000, 33000000, 35000000, 29000000)
#' )
#' hla_data <- HLA_get(example_data, chromosome_col="CHR", position_col="POS")
#' print(hla_data)
HLA_get <- function(data, chromosome_col="CHR", position_col="BP", lower_bound=25e6, upper_bound=34e6) {
  hla_data <- data %>%
    dplyr::filter(!!sym(chromosome_col) == 6 &
                    !!sym(position_col) >= lower_bound &
                    !!sym(position_col) <= upper_bound)
  return(hla_data)
}
