# conditional analysis ----
#' locate the significant SNP for conditional analysis
#'
#' @param x data.frame of the SNP information
#' @param environment.list tibble of the environment list (see the example for usage)
#' @param significance_level numeric, significance level, default is 5e-8
#'
#' @return message that informs the user of the independant SNP for the following conditional analysis
#' @export
#' @section HLA data analysis:
#'
#' @examples
#' con_dir <- "/Users/leoarrow/project/VKH2024/data/zuo/con_su" # specify the directory to store the HLA original data and subsequent conditional analysis results data
#' files <- list.files(con_dir,full.names = T) %>% as.vector(); files # update it each time
#' x1 <- fread(files[1]) %>% arrange(desc(CHISQ)) # read the data and sort it by CHISQ/P value
#' x2 <- fread(files[2]) %>% arrange(P) # repeat it until no more independent signal can be found.
#' x3 <- fread(files[3]) %>% arrange(P)
#' x4 <- fread(files[4]) %>% arrange(P)
#' x5 <- fread(files[5]) %>% arrange(P)
#' x6 <- fread(files[6]) %>% arrange(P)
#' x7 <- fread(files[7]) %>% arrange(P)
#' x8 <- fread(files[8]) %>% arrange(P)
#'
#' env <- ls() # get the environment; this line if were put in main func will lead to error.
#' environment.list <- tibble(item = as.vector(grep("^x[0-9]+$", x = env, value = TRUE))) # load the environment
#' check_significant_SNP(x8, environment.list, significance_level = 5e-8)
#'
check_significant_SNP <- function(x, environment.list, significance_level = 5e-8) {
  if (nrow(environment.list) == 0) {
    message("No matching variables found in the environment.")
    return()
  }
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
      break
    } else {
      leo_message("Not this one! Let's move on to the next one.")
    }

    if (i == nrow(x)) {
      leo_message("The conditional analysis can be over now. So long!")
    }
  }
}

