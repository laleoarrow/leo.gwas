#' Title Get ID using CHR, BP, A2 (REF/Non-effect), A1 (ALT/Effect)
#'
#' @param x df data.frame with CHR, BP, A2, A1 columns
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