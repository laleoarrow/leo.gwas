.onLoad <- function(libname, pkgname) {
  ascii_arts <- c(
    "
  ==============================================
            \u2588\u2588\u2557     \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2557
            \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255d\u2588\u2588\u2554\u2550\u2550\u2550\u2588\u2588\u2557
            \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2551   \u2588\u2588\u2551
            \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u255d  \u2588\u2588\u2551   \u2588\u2588\u2551
            \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u255a\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255d
            \u255a\u2550\u2550\u2550\u2550\u2550\u2550\u255d\u255a\u2550\u2550\u2550\u2550\u2550\u2550\u255d \u255a\u2550\u2550\u2550\u2550\u2550\u255d
  ==============================================
  ",
    "
  =============================================
            ___       ___           ___
           /\\__\\     /\\  \\         /\\  \\
          /:/  /    /::\\  \\       /::\\  \\
         /:/  /    /:/\\:\\  \\     /:/\\:\\  \\
        /:/  /    /::\\~\\:\\  \\   /:/  \\:\\  \\
       /:/__/    /:/\\:\\ \\:\\__\\ /:/__/ \\:\\__\\
       \\:\\  \\    \\:\\~\\:\\ \\/__/ \\:\\  \\ /:/  /
        \\:\\  \\    \\:\\ \\:\\__\\    \\:\\  /:/  /
         \\:\\  \\    \\:\\ \\/__/     \\:\\/:/  /
          \\:\\__\\    \\:\\__\\        \\::/  /
           \\/__/     \\/__/         \\/__/
  =============================================
    "
  )

  .leo_startup_message(sample(ascii_arts, 1))
  .leo_startup_message(" \U0001F981\U0001F981\U0001F981 Welcome to use the LEO package ! \U0001F981\U0001F981\U0001F981")
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .leo_startup_message(paste0(" \U0001F981\U0001F981\U0001F981 System Time: ", current_time, " \U0001F981\U0001F981\U0001F981"))
}

.leo_startup_message <- function(text) {
  if (.leo_supports_truecolor()) {
    message(.leo_gradient(text))
  } else {
    leo_message(text, color = "33")
  }
}

.leo_supports_truecolor <- function() {
  if (nzchar(Sys.getenv("NO_COLOR"))) {
    return(FALSE)
  }
  if (nzchar(Sys.getenv("RSTUDIO"))) {
    return(FALSE)
  }
  term <- Sys.getenv("TERM")
  if (identical(term, "dumb")) {
    return(FALSE)
  }
  colorterm <- Sys.getenv("COLORTERM")
  if (nzchar(colorterm) && grepl("truecolor|24bit", colorterm, ignore.case = TRUE)) {
    return(TRUE)
  }
  FALSE
}

.leo_gradient <- function(text,
                          start = c(255, 240, 120),
                          end = c(255, 180, 0)) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  colored <- vapply(lines, function(line) {
    if (!nzchar(line)) {
      return("")
    }

    chars <- strsplit(line, "", fixed = TRUE)[[1]]
    n <- length(chars)
    if (n == 1L) {
      t <- 0
    } else {
      t <- (seq_len(n) - 1) / (n - 1)
    }

    r <- round(start[1] + (end[1] - start[1]) * t)
    g <- round(start[2] + (end[2] - start[2]) * t)
    b <- round(start[3] + (end[3] - start[3]) * t)

    paste0(
      paste0(sprintf("\033[38;2;%d;%d;%dm%s", r, g, b, chars), collapse = ""),
      "\033[0m"
    )
  }, character(1))

  paste(colored, collapse = "\n")
}
