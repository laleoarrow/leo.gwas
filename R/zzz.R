.onLoad <- function(libname, pkgname) {
  leo_message(" 🦁🦁🦁 Welcome to use the LEO package ! 🦁🦁🦁", color = "33")
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  leo_message(paste0(" 🦁🦁🦁 System Time: ", current_time, " 🦁🦁🦁"), color = "33")
}
