# Internal helper: run a shell command and stop on error.
.csmr_run <- function(cmd, args = character(), env = character()) {
  leo.basic::leo_log("Run: {cmd} {paste(shQuote(args), collapse = ' ')}")
  out <- suppressWarnings(system2(cmd, args = args, env = env, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (!identical(status, 0L)) {
    tail_out <- utils::tail(out, 20L)
    stop(paste0("Command failed: ", cmd, " ", paste(shQuote(args), collapse = " "),
                "\n", paste(tail_out, collapse = "\n")))
  }
  out
}

# Internal helper: get a valid conda executable path.
.csmr_conda <- function(conda = Sys.which("conda")) {
  if (!nzchar(conda)) stop("`conda` was not found. Please provide a valid conda path or add it to PATH.")
  normalizePath(conda, winslash = "/", mustWork = TRUE)
}

# Internal helper: get env prefix path from `conda env list`.
.csmr_env_prefix <- function(conda, env_name) {
  out <- .csmr_run(conda, c("env", "list"))
  x <- trimws(out)
  x <- x[nzchar(x) & !grepl("^#|^\\*", x)]
  hit <- x[grepl(paste0("^", env_name, "\\s"), x)]
  if (length(hit) == 0L) return("")
  tokens <- strsplit(hit[1], "\\s+")[[1]]
  as.character(utils::tail(tokens, 1))
}

# Internal helper: write a temporary compatibility env file.
.csmr_compat_env_file <- function(env_name) {
  deps <- c(
    "python=3.10",
    "snakemake=7.18.2",
    "bedtools=2.31.1",
    "pybedtools=0.12.0",
    "r-base=4.5",
    "r-getopt",
    "r-coloc",
    "r-susier",
    "r-pkgload",
    "r-ieugwasr",
    "r-twosamplemr",
    "r-data.table",
    "r-dplyr",
    "r-remotes",
    "plink=1.90"
  )
  yml <- c(
    paste0("name: ", env_name),
    "channels:",
    "  - conda-forge",
    "  - bioconda",
    "dependencies:",
    paste0("  - ", deps)
  )
  tf <- tempfile(fileext = ".yml")
  writeLines(yml, tf)
  tf
}

# Internal helper: install the extra R packages used by csMR scripts.
.csmr_install_r_deps <- function(conda, env_name, r_lib) {
  script <- c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "r_lib <- normalizePath(args[[1]], winslash = '/', mustWork = FALSE)",
    "dir.create(r_lib, recursive = TRUE, showWarnings = FALSE)",
    "Sys.setenv(R_LIBS = r_lib)",
    ".libPaths(r_lib)",
    "options(timeout = max(1200, getOption('timeout')))",
    "options(repos = c(CRAN = 'https://cloud.r-project.org'))",
    "locks <- Sys.glob(file.path(r_lib, '00LOCK*'))",
    "if (length(locks) > 0L) unlink(locks, recursive = TRUE, force = TRUE)",
    "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', lib = r_lib)",
    "install_gh <- function(pkg, repo, ref = NULL, deps = TRUE) {",
    "  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))",
    "  args <- list(repo = repo, lib = r_lib, upgrade = 'never', force = TRUE, dependencies = deps)",
    "  if (!is.null(ref)) args$ref <- ref",
    "  do.call(remotes::install_github, args)",
    "  invisible(TRUE)",
    "}",
    "install_gh('phenoscanner', 'phenoscanner/phenoscanner')",
    "if (!requireNamespace('mr.raps', quietly = TRUE)) {",
    "  fit <- try(install_gh('mr.raps', 'qingyuanzhao/mr.raps', ref = '0.4.1'), silent = TRUE)",
    "  if (inherits(fit, 'try-error')) install_gh('mr.raps', 'qingyuanzhao/mr.raps', ref = 'v0.4.1')",
    "}",
    "install_gh('RadialMR', 'WSpiller/RadialMR')",
    "install_gh('MRPRESSO', 'rondolab/MR-PRESSO')",
    "if (!requireNamespace('TwoSampleMR', quietly = TRUE)) install_gh('TwoSampleMR', 'MRCIEU/TwoSampleMR', deps = FALSE)",
    "need <- c('data.table', 'dplyr', 'getopt', 'coloc', 'ieugwasr', 'phenoscanner', 'mr.raps', 'RadialMR', 'MRPRESSO', 'TwoSampleMR')",
    "miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]",
    "if (length(miss) > 0L) stop('Missing R packages: ', paste(miss, collapse = ', '))",
    "cat('R_LIBS=', Sys.getenv('R_LIBS'), '\\n', sep = '')",
    "print(.libPaths())"
  )
  tf <- tempfile(fileext = ".R")
  on.exit(unlink(tf), add = TRUE)
  writeLines(script, tf)
  .csmr_run(conda, c("run", "--no-capture-output", "-n", env_name, "Rscript", tf, r_lib), env = paste0("R_LIBS=", r_lib))
}

# Internal helper: verify the final R package set in a fresh conda R session.
.csmr_check_r_deps <- function(conda, env_name, r_lib) {
  script <- c(
    "pkgs <- c('data.table', 'dplyr', 'getopt', 'coloc', 'ieugwasr', 'phenoscanner', 'mr.raps', 'RadialMR', 'MRPRESSO', 'TwoSampleMR')",
    "status <- setNames(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE), pkgs)",
    "print(status)",
    "print(.libPaths())",
    "cat(Sys.getenv('R_LIBS'), '\\n')",
    "if (!all(status)) stop('Missing R packages after setup: ', paste(names(status)[!status], collapse = ', '))"
  )
  tf <- tempfile(fileext = ".R")
  on.exit(unlink(tf), add = TRUE)
  writeLines(script, tf)
  .csmr_run(conda, c("run", "--no-capture-output", "-n", env_name, "Rscript", tf), env = paste0("R_LIBS=", r_lib))
}

# Internal helper: normalize supported data input.
.csmr_read_input <- function(data) {
  if (is.data.frame(data)) return(as.data.frame(data, stringsAsFactors = FALSE))
  if (!is.character(data) || length(data) != 1L) stop("`data` must be a data.frame or a file path.")
  data.table::fread(path.expand(data), data.table = FALSE, showProgress = FALSE)
}

# Internal helper: make a stable id from file name.
.csmr_make_id <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.(tsv|txt|ma|csv)(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  gsub("--+", "-", x)
}

# Internal helper: convert one table to csMR format.
.csmr_format_df <- function(df, type, SNP_col, GENE_col, A1_col, A2_col, MAF_col, BETA_col, SE_col, P_col, N_col, n = NULL) {
  type <- match.arg(type, c("gwas", "qtl"))
  need <- c(SNP_col, A1_col, A2_col, MAF_col, BETA_col, SE_col, P_col, N_col)
  out_names <- c("SNP", "A1", "A2", "MAF", "BETA", "SE", "P", "N")
  if (identical(type, "qtl")) {
    need <- c(SNP_col, GENE_col, A1_col, A2_col, MAF_col, BETA_col, SE_col, P_col, N_col)
    out_names <- c("SNP", "GENE", "A1", "A2", "MAF", "BETA", "SE", "P", "N")
  }
  if (!all(need %in% names(df))) stop("Missing required columns in `data`.")
  df <- df[, need, drop = FALSE]
  names(df) <- out_names
  if (identical(type, "qtl")) df$GENE <- as.character(df$GENE)
  df$SNP <- as.character(df$SNP)
  df$A1 <- toupper(as.character(df$A1))
  df$A2 <- toupper(as.character(df$A2))
  df <- get_biallelic_snp(df)
  leo.basic::leo_log("keep only biallelic SNPs.")
  df$MAF <- suppressWarnings(as.numeric(df$MAF))
  df$BETA <- suppressWarnings(as.numeric(df$BETA))
  df$SE <- suppressWarnings(as.numeric(df$SE))
  df$P <- suppressWarnings(as.numeric(df$P))
  df$N <- suppressWarnings(as.numeric(df$N))
  if (identical(type, "gwas")) {
    miss_n <- !is.finite(df$N) | df$N <= 0
    if (any(miss_n)) {
      if (is.null(n)) stop("Missing `N` values detected; provide `n` to fill them.")
      df$N[miss_n] <- as.numeric(n)
    }
  }
  df$MAF <- pmin(df$MAF, 1 - df$MAF)
  keep <- nzchar(df$SNP) & nzchar(df$A1) & nzchar(df$A2)
  if (identical(type, "qtl")) keep <- keep & nzchar(df$GENE)
  keep <- keep & grepl("^[ACGT]+$", df$A1) & grepl("^[ACGT]+$", df$A2)
  keep <- keep & is.finite(df$MAF) & df$MAF >= 0.01 & df$MAF < 1
  keep <- keep & is.finite(df$BETA) & is.finite(df$SE) & df$SE > 0
  keep <- keep & is.finite(df$P) & df$P > 0 & df$P < 1
  keep <- keep & is.finite(df$N) & df$N > 0
  df[keep, , drop = FALSE]
}

# Internal helper: resolve qtl input for step2.
.csmr_qtl_paths <- function(qtl_input) {
  if (is.data.frame(qtl_input)) {
    path_col <- intersect(c("output", "path"), names(qtl_input))
    if (length(path_col) == 0L) stop("`qtl_input` data.frame must contain `output` or `path`.")
    return(normalizePath(path.expand(as.character(qtl_input[[path_col[1]]])), winslash = "/", mustWork = TRUE))
  }
  if (is.character(qtl_input) && length(qtl_input) == 1L && dir.exists(path.expand(qtl_input))) {
    hit <- list.files(path.expand(qtl_input), pattern = "\\.ma$", full.names = TRUE)
    if (length(hit) == 0L) stop("No `.ma` files found in `qtl_input` directory.")
    return(normalizePath(hit, winslash = "/", mustWork = TRUE))
  }
  if (is.character(qtl_input) && length(qtl_input) > 0L) {
    return(normalizePath(path.expand(qtl_input), winslash = "/", mustWork = TRUE))
  }
  stop("`qtl_input` must be a manifest data.frame, a `.ma` directory, or a `.ma` path vector.")
}

# Internal helper: parse base output dir and ref path from config.yml.
.csmr_parse_config <- function(config_file) {
  lines <- readLines(config_file, warn = FALSE)
  base_line <- grep("^\\s*BASE_OUTPUT_DIR\\s*:", lines, value = TRUE)
  base_output_dir <- trimws(sub("^\\s*BASE_OUTPUT_DIR\\s*:\\s*", "", base_line[1]))
  base_output_dir <- gsub("^['\"]|['\"]$", "", base_output_dir)
  idx <- grep("^\\s*GWAS_REFERENCE_GENOTYPE\\s*:", lines)
  ref_path <- ""
  if (length(idx) > 0L) {
    ref_line <- grep("^\\s*path\\s*:", lines[(idx[1] + 1L):length(lines)], value = TRUE)
    if (length(ref_line) > 0L) ref_path <- trimws(sub("^\\s*path\\s*:\\s*", "", ref_line[1]))
  }
  ref_path <- gsub("^['\"]|['\"]$", "", ref_path)
  outcome_line <- grep("^\\s*OUTCOME_DIR\\s*:", lines, value = TRUE)
  outcome_dir <- ""
  if (length(outcome_line) > 0L) {
    outcome_dir <- trimws(sub("^\\s*OUTCOME_DIR\\s*:\\s*", "", outcome_line[1]))
    outcome_dir <- gsub("^['\"]|['\"]$", "", outcome_dir)
  }
  list(base_output_dir = normalizePath(path.expand(base_output_dir), winslash = "/", mustWork = FALSE),
       ref_genotype_dir = normalizePath(path.expand(ref_path), winslash = "/", mustWork = TRUE),
       outcome_dir = if (nzchar(outcome_dir)) normalizePath(path.expand(outcome_dir), winslash = "/", mustWork = FALSE) else "")
}

#' Configure Conda Environment for csMR
#'
#' @description
#' `csMR_env()` prepares the official csMR environment. It clones the official
#' repository if missing, tries `envs/envpy3.yml` first, falls back to a
#' temporary compatibility env file only when the official file cannot solve,
#' installs required R packages into the csMR conda library, checks PLINK 1.90,
#' and returns the runtime paths required by the csMR README.
#'
#' @param repo_dir Path to the csMR repository.
#' @param repo_url Official csMR GitHub URL.
#' @param env_name Conda environment name.
#' @param conda Path to conda.
#' @param mamba Path to mamba.
#' @param plink_path Optional PLINK 1.90 path. If empty, `PATH` and
#'   `plinkbinr::get_plink_exe()` are tried.
#' @param ref Git branch or tag.
#' @param overwrite Whether to rebuild an existing env.
#' @return A list with `repo_dir`, `env_name`, `env_prefix`, `env_file`,
#'   `r_home`, `r_library`, and `plink`.
#' @examples
#' \dontrun{
#' cfg <- csMR_env()
#' cfg$r_library
#' }
#' @importFrom cli cat_rule
#' @importFrom cli cli_alert_info
#' @export
csMR_env <- function(
    repo_dir = "~/Project/software/csMR",
    repo_url = "https://github.com/rhhao/csMR.git",
    env_name = "csMR",
    conda = Sys.which("conda"),
    mamba = Sys.which("mamba"),
    plink_path = Sys.which("plink"),
    ref = "main",
    overwrite = FALSE
) {
  cli::cat_rule("Preparing csMR environment ...", line_col = "blue", col = "blue")
  conda <- .csmr_conda(conda)
  repo_dir <- path.expand(repo_dir)
  if (basename(repo_dir) != "csMR") repo_dir <- file.path(repo_dir, "csMR")
  solver <- conda
  if (nzchar(mamba)) solver <- normalizePath(mamba, winslash = "/", mustWork = TRUE)
  if (!dir.exists(repo_dir)) {
    git_bin <- Sys.which("git")
    if (!nzchar(git_bin)) stop("`git` was not found in PATH.")
    dir.create(dirname(repo_dir), recursive = TRUE, showWarnings = FALSE)
    .csmr_run(git_bin, c("clone", "--depth", "1", "--branch", ref, repo_url, repo_dir))
  }
  official_env_file <- file.path(repo_dir, "envs", "envpy3.yml")
  if (!file.exists(official_env_file)) stop("Cannot find csMR env file: ", official_env_file)
  env_prefix <- .csmr_env_prefix(conda, env_name)
  env_file <- official_env_file
  if (nzchar(env_prefix) && isTRUE(overwrite)) {
    .csmr_run(conda, c("env", "remove", "-n", env_name, "--yes"))
    env_prefix <- ""
  }
  if (!nzchar(env_prefix)) {
    cli::cli_alert_info("Creating conda env from official csMR env file ...")
    fit <- try(.csmr_run(solver, c("env", "create", "--name", env_name, "--file", official_env_file, "--yes")), silent = TRUE)
    if (inherits(fit, "try-error")) {
      leo.basic::leo_log("Official envpy3.yml failed on this machine; retry with a temporary compatibility env file.", level = "warning")
      if (nzchar(.csmr_env_prefix(conda, env_name))) .csmr_run(conda, c("env", "remove", "-n", env_name, "--yes"))
      env_file <- .csmr_compat_env_file(env_name)
      on.exit(unlink(env_file), add = TRUE)
      .csmr_run(solver, c("env", "create", "--name", env_name, "--file", env_file, "--yes"))
    }
    env_prefix <- .csmr_env_prefix(conda, env_name)
  } else {
    leo.basic::leo_log("Conda env `{env_name}` already exists; skip creation.", level = "info")
  }
  if (!nzchar(env_prefix)) stop("Failed to determine conda env prefix for: ", env_name)
  r_lib <- file.path(env_prefix, "lib", "R", "library")
  r_home <- file.path(env_prefix, "lib", "R")
  r_env <- paste0("R_LIBS=", r_lib)

  cli::cli_alert_info("Checking PLINK availability ...")
  has_plink <- suppressWarnings(system2(conda, c("run", "-n", env_name, "plink", "--version"), stdout = TRUE, stderr = TRUE))
  if (!any(grepl("PLINK v1\\.9(0|\\b)", has_plink))) {
    plink_bin <- ""
    if (nzchar(plink_path)) plink_bin <- normalizePath(path.expand(plink_path), winslash = "/", mustWork = TRUE)
    if (!nzchar(plink_bin)) plink_bin <- Sys.which("plink")
    if (!nzchar(plink_bin) && requireNamespace("plinkbinr", quietly = TRUE)) {
      plink_bin <- tryCatch(plinkbinr::get_plink_exe(), error = function(e) "")
    }
    if (!nzchar(plink_bin)) stop("`plink` was not found in PATH or `plinkbinr`.")
    plink_bin <- normalizePath(path.expand(plink_bin), winslash = "/", mustWork = TRUE)
    plink_version <- paste(suppressWarnings(system2(plink_bin, "--version", stdout = TRUE, stderr = TRUE)), collapse = "\n")
    if (!grepl("PLINK v1\\.9(0|\\b)", plink_version)) stop("`plink` must be version 1.9x.")
    dir.create(file.path(env_prefix, "bin"), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(plink_bin, file.path(env_prefix, "bin", "plink"), overwrite = TRUE)
    if (!isTRUE(ok)) stop("Failed to copy `plink` into the csMR conda env.")
    Sys.chmod(file.path(env_prefix, "bin", "plink"), mode = "0755")
  }

  cli::cli_alert_info("Installing required R packages for csMR ...")
  .csmr_install_r_deps(conda, env_name, r_lib)
  cli::cli_alert_info("Running final environment checks ...")
  .csmr_check_r_deps(conda, env_name, r_lib)
  .csmr_run(conda, c("run", "-n", env_name, "snakemake", "--version"), env = r_env)
  .csmr_run(conda, c("run", "-n", env_name, "plink", "--version"), env = r_env)
  plink_exec <- file.path(env_prefix, "bin", "plink")
  if (!file.exists(plink_exec)) {
    plink_exec <- suppressWarnings(system2(conda, c("run", "-n", env_name, "which", "plink"), stdout = TRUE, stderr = TRUE))
    status <- attr(plink_exec, "status")
    if (is.null(status)) status <- 0L
    if (!identical(status, 0L)) plink_exec <- ""
    plink_exec <- trimws(plink_exec[1])
  }
  if (!nzchar(plink_exec) || !file.exists(plink_exec)) stop("Failed to resolve `plink` path in `csMR` env.")

  out <- list(
    repo_dir = normalizePath(repo_dir, winslash = "/", mustWork = TRUE),
    env_name = env_name,
    env_prefix = normalizePath(env_prefix, winslash = "/", mustWork = TRUE),
    env_file = env_file,
    r_home = normalizePath(r_home, winslash = "/", mustWork = FALSE),
    r_library = normalizePath(r_lib, winslash = "/", mustWork = FALSE),
    plink = normalizePath(plink_exec, winslash = "/", mustWork = TRUE)
  )
  leo.basic::leo_log("csMR environment setup completed.", level = "success")
  out
}

#' Prepare csMR Step1 Input
#'
#' @description
#' Convert GWAS or QTL summary statistics to the csMR-required `.ma` format
#' using explicit column mappings. This function only performs thin formatting
#' and basic QC; it does not guess genome build or map non-rsID SNPs.
#'
#' @param data GWAS/QTL data.frame, file path, file vector, or a QTL directory.
#' @param type `"gwas"` or `"qtl"`.
#' @param output Output `.ma` file path for GWAS or single-table QTL, or output
#'   directory for multi-file QTL input.
#' @param SNP_col Input SNP column name.
#' @param GENE_col Input gene column name for QTL.
#' @param A1_col Input effect allele column name.
#' @param A2_col Input other allele column name.
#' @param MAF_col Input MAF column name.
#' @param BETA_col Input beta column name.
#' @param SE_col Input SE column name.
#' @param P_col Input P column name.
#' @param N_col Input sample size column name.
#' @param n Optional fixed sample size used to fill missing `N` in GWAS.
#' @return For GWAS, a list with `output`, `n_input`, `n_output`, and
#'   `n_missing_filled`. For QTL, a manifest data.frame with `id`, `input`,
#'   `output`, `n_input`, and `n_output`.
#' @examples
#' \dontrun{
#' csMR_step1_prep(
#'   data = "~/Project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv",
#'   type = "gwas",
#'   output = "~/Project/iridocyclitis/output/csMR/step1/exposure.ma",
#'   SNP_col = "rsID",
#'   A1_col = "effect_allele",
#'   A2_col = "other_allele",
#'   MAF_col = "EAF",
#'   BETA_col = "beta",
#'   SE_col = "se",
#'   P_col = "pval",
#'   N_col = "N"
#' )
#' }
#' @importFrom cli cat_rule
#' @importFrom cli cli_alert_info
#' @export
csMR_step1_prep <- function(
    data,
    type = c("gwas", "qtl"),
    output,
    SNP_col = "SNP",
    GENE_col = "GENE",
    A1_col = "A1",
    A2_col = "A2",
    MAF_col = "MAF",
    BETA_col = "BETA",
    SE_col = "SE",
    P_col = "P",
    N_col = "N",
    n = NULL
) {
  type <- match.arg(type)
  cli::cat_rule("Preparing csMR step1 inputs ...", line_col = "blue", col = "blue")
  cli::cli_alert_info("Formatting `{type}` input for csMR.")
  if (identical(type, "gwas")) {
    df <- .csmr_read_input(data)
    n_input <- nrow(df)
    miss_n <- sum(!is.finite(suppressWarnings(as.numeric(df[[N_col]]))) | suppressWarnings(as.numeric(df[[N_col]])) <= 0, na.rm = TRUE)
    df <- .csmr_format_df(df, type = "gwas", SNP_col = SNP_col, GENE_col = GENE_col, A1_col = A1_col,
                          A2_col = A2_col, MAF_col = MAF_col, BETA_col = BETA_col, SE_col = SE_col,
                          P_col = P_col, N_col = N_col, n = n)
    output <- path.expand(output)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
    out <- list(output = normalizePath(output, winslash = "/", mustWork = FALSE), n_input = n_input,
                n_output = nrow(df), n_missing_filled = as.integer(miss_n))
    leo.basic::leo_log("Prepared csMR GWAS file with {out$n_output} row{?s}.", level = "success")
    return(out)
  }

  qtl_files <- data
  if (is.character(qtl_files) && length(qtl_files) == 1L && dir.exists(path.expand(qtl_files))) {
    qtl_files <- list.files(path.expand(qtl_files), pattern = "\\.(tsv|txt|csv|ma)(\\.gz)?$", full.names = TRUE)
  }
  if (is.data.frame(qtl_files)) {
    df <- .csmr_format_df(qtl_files, type = "qtl", SNP_col = SNP_col, GENE_col = GENE_col, A1_col = A1_col,
                          A2_col = A2_col, MAF_col = MAF_col, BETA_col = BETA_col, SE_col = SE_col,
                          P_col = P_col, N_col = N_col)
    output <- path.expand(output)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
    out <- data.frame(id = .csmr_make_id(output), input = "data.frame",
                      output = normalizePath(output, winslash = "/", mustWork = FALSE),
                      n_input = nrow(qtl_files), n_output = nrow(df), stringsAsFactors = FALSE)
    leo.basic::leo_log("Prepared 1 csMR QTL file.", level = "success")
    return(out)
  }
  if (!is.character(qtl_files) || length(qtl_files) == 0L) stop("`data` must be a QTL file path, file vector, directory, or data.frame.")
  qtl_files <- normalizePath(path.expand(qtl_files), winslash = "/", mustWork = TRUE)
  output <- path.expand(output)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  out <- vector("list", length(qtl_files))
  for (i in seq_along(qtl_files)) {
    df_raw <- .csmr_read_input(qtl_files[i])
    df <- .csmr_format_df(df_raw, type = "qtl", SNP_col = SNP_col, GENE_col = GENE_col, A1_col = A1_col,
                          A2_col = A2_col, MAF_col = MAF_col, BETA_col = BETA_col, SE_col = SE_col,
                          P_col = P_col, N_col = N_col)
    output_file <- file.path(output, paste0(.csmr_make_id(qtl_files[i]), ".ma"))
    data.table::fwrite(df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    out[[i]] <- data.frame(id = .csmr_make_id(qtl_files[i]), input = qtl_files[i],
                           output = normalizePath(output_file, winslash = "/", mustWork = FALSE),
                           n_input = nrow(df_raw), n_output = nrow(df), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, out)
  leo.basic::leo_log("Prepared {nrow(out)} csMR QTL file{?s}.", level = "success")
  out
}

#' Build csMR config.yml
#'
#' @description
#' Write an official-style csMR `config.yml` for Step 3.
#'
#' @details
#' Uppercase arguments such as `BASE_OUTPUT_DIR`, `GWAS_REFERENCE_GENOTYPE`,
#' and `eQTL_REFERENCE_GENOTYPE` intentionally mirror official csMR
#' `config.yml` keys to keep mapping explicit and reduce confusion when
#' debugging Step 3.
#'
#' @param repo_dir csMR repository directory.
#' @param config.yml_to Output config file path.
#' @param BASE_OUTPUT_DIR Output root used by csMR Step 3.
#' @param qtl_input_dir Directory of QTL `.ma` files (or any input supported by
#'   `.csmr_qtl_paths()`).
#' @param exposure_ma Exposure GWAS `.ma` file path.
#' @param exposure_id Exposure id in `GWAS_SUMSTATS`.
#' @param exposure_type Exposure type, either `"cc"` or `"quant"`.
#' @param exposure_ma_dir Optional directory of multiple exposure `.ma` files.
#'   If provided, `exposure_ma` is ignored and ids are generated from file names.
#' @param outcome_ma_dir Outcome `.ma` directory (used as `OUTCOME_DIR`).
#'   Must contain **only** `.ma` files — Snakemake's `os.listdir()` picks up
#'   every item (incl. `.DS_Store`, subdirectories) as an outcome wildcard.
#' @param GWAS_REFERENCE_GENOTYPE GWAS reference genotype directory. If `NULL`,
#'   default to `repo_dir/data/reference_genome_1000G_EUR`.
#' @param eQTL_REFERENCE_GENOTYPE eQTL reference genotype directory. If `NULL`,
#'   default to `repo_dir/data/reference_genome_1000G_EUR`.
#' @param duplicated_snp_path Path to duplicated SNP file for the reference
#'   genotype. Use `"None"` if not applicable (official csMR default).
#' @param coloc_window_size_bp Coloc window size.
#' @param coloc_coverages Coloc coverage.
#' @param coloc_threads Coloc threads.
#' @param coloc_cutoff Coloc cutoff.
#' @return Absolute path to the written config file.
#' @examples
#' \dontrun{
#' out_cfg <- csMR_step2_config.yml(
#'   repo_dir = "~/Project/software/csMR",
#'   config.yml_to = "./output/csMR/step2_config/config.yml",
#'   BASE_OUTPUT_DIR = "./output/csMR/step3_run", # where to store the final csMR output
#'   qtl_input_dir = "./output/csMR/step1_data_preparation/sc_qtl_dir1",
#'   exposure_ma = "./output/csMR/step1_data_preparation/gwas_exp_p1.ma",
#'   exposure_id = "exp_p1",
#'   exposure_type = "cc",
#'   outcome_ma_dir = "./output/csMR/step1_data_preparation/outcome"
#' )
#' # Example console messages:
#' # i [23:36:24] Writing csMR config.yml ...
#' # v [23:36:24] csMR step2 config written: ./output/csMR/step2_config/config.yml
#' out_cfg
#' # [1] "/Users/leoarrow/Project/iridocyclitis/output/csMR/step2_config/config.yml"
#' }
#' @export
csMR_step2_config.yml <- function(
    repo_dir = "~/Project/software/csMR",
    config.yml_to = "./output/csMR/step2_config/config.yml",
    BASE_OUTPUT_DIR = "./output/csMR/step3_run", # where to store the final csMR output
    qtl_input_dir = "./output/csMR/step1_data_preparation/sc_qtl_dir1",
    exposure_ma = "./output/csMR/step1_data_preparation/gwas_exp_p1.ma",
    exposure_id = "exp_p1",
    exposure_type = "cc",
    exposure_ma_dir = NULL,
    outcome_ma_dir = "./output/csMR/step1_data_preparation/outcome",
    GWAS_REFERENCE_GENOTYPE = NULL,
    eQTL_REFERENCE_GENOTYPE = NULL,
    duplicated_snp_path = "None",
    coloc_window_size_bp = 100000L,
    coloc_coverages = 0.9,
    coloc_threads = 5L,
    coloc_cutoff = 0.8
) {
  leo.basic::leo_log("Writing csMR config.yml ...", level = "info")
  repo_dir <- normalizePath(path.expand(repo_dir), winslash = "/", mustWork = TRUE)
  if (!exposure_type %in% c("cc", "quant")) stop("`exposure_type` must be `cc` or `quant`.")
  if (is.null(GWAS_REFERENCE_GENOTYPE)) GWAS_REFERENCE_GENOTYPE <- file.path(repo_dir, "data", "reference_genome_1000G_EUR")
  if (is.null(eQTL_REFERENCE_GENOTYPE)) eQTL_REFERENCE_GENOTYPE <- file.path(repo_dir, "data", "reference_genome_1000G_EUR")

  qtl_paths <- .csmr_qtl_paths(qtl_input_dir)

  if (!is.null(exposure_ma_dir)) {
    exposure_paths <- list.files(path.expand(exposure_ma_dir), pattern = "\\.ma$", full.names = TRUE)
    if (length(exposure_paths) == 0L) stop("No `.ma` files found in `exposure_ma_dir`.")
  } else {
    if (is.list(exposure_ma) && "output" %in% names(exposure_ma)) exposure_ma <- exposure_ma$output
    exposure_paths <- exposure_ma
  }
  exposure_paths <- normalizePath(path.expand(exposure_paths), winslash = "/", mustWork = TRUE)

  outcome_dir <- normalizePath(path.expand(outcome_ma_dir), winslash = "/", mustWork = TRUE)
  # Auto-remove .DS_Store (macOS Finder artifact) — always safe, prevents Snakemake from treating it as an outcome
  ds_store <- file.path(outcome_dir, ".DS_Store")
  if (file.exists(ds_store)) {
    unlink(ds_store)
    leo.basic::leo_log("Removed .DS_Store from outcome_ma_dir.", level = "info")
  }
  outcome_files <- list.files(outcome_dir, pattern = "\\.ma$", full.names = TRUE)
  if (length(outcome_files) == 0L) stop("No `.ma` files found in `outcome_ma_dir`.")
  all_items <- list.files(outcome_dir, all.files = FALSE)
  non_ma <- setdiff(all_items, list.files(outcome_dir, pattern = "\\.ma$"))
  if (length(non_ma) > 0L) {
    stop(
      "outcome_ma_dir contains non-.ma items: ", paste(non_ma, collapse = ", "),
      ". Snakemake's os.listdir() will treat ALL items as outcomes. ",
      "Please use a directory containing ONLY .ma files."
    )
  }
  qtl_ids <- vapply(qtl_paths, .csmr_make_id, character(1))
  all_ids <- if (!is.null(exposure_ma_dir)) {
    c(qtl_ids, vapply(exposure_paths, .csmr_make_id, character(1)))
  } else {
    c(qtl_ids, exposure_id)
  }
  bad_ids <- all_ids[grepl("--", all_ids, fixed = TRUE)]
  if (length(bad_ids) > 0L) stop("IDs must not contain double dashes (--): ", paste(bad_ids, collapse = ", "))

  gwas_ref <- normalizePath(path.expand(GWAS_REFERENCE_GENOTYPE), winslash = "/", mustWork = TRUE)
  eqtl_ref <- normalizePath(path.expand(eQTL_REFERENCE_GENOTYPE), winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(path.expand(BASE_OUTPUT_DIR), winslash = "/", mustWork = FALSE)
  config_file <- path.expand(config.yml_to)
  dir.create(dirname(config_file), recursive = TRUE, showWarnings = FALSE)
  yaml_quote <- function(x) paste0("\"", gsub("\"", "\\\\\"", as.character(x)), "\"")
  lines <- c(
    "---",
    "",
    "# All file paths use absolute paths. Generated by leo.gwas::csMR_step2_config.yml()",
    "# Environment variables (including ~) are not processed by csMR.",
    "",
    "# Folder where all the outputs are saved (COLOC, MR results files, etc.)",
    paste0("BASE_OUTPUT_DIR: ", yaml_quote(output_dir)),
    "",
    "### eQTL_INPUT: list of eQTL files",
    "# id: identifier for the dataset. *MUST be unique and NO double dash allowed.",
    "# path: filepath to the QTL data. *Compressed files are not allowed.",
    "eQTL_INPUT:")
  for (i in seq_along(qtl_paths)) {
    lines <- c(lines,
               paste0("  - id: ", yaml_quote(.csmr_make_id(qtl_paths[i]))),
               paste0("    path: ", yaml_quote(qtl_paths[i])))
  }
  lines <- c(lines,
             "",
             "### GWAS_SUMSTATS:",
             "# id: identifier for GWAS. *MUST be unique and NO double dash allowed.",
             "# path: filepath to GWAS sumstats. *Compressed files are not allowed.",
             "GWAS_SUMSTATS:")
  if (!is.null(exposure_ma_dir)) {
    for (i in seq_along(exposure_paths)) {
      lines <- c(lines,
                 paste0("  - id: ", yaml_quote(.csmr_make_id(exposure_paths[i]))),
                 paste0("    path: ", yaml_quote(exposure_paths[i])),
                 paste0("    type: ", exposure_type))
    }
  } else {
    lines <- c(lines,
               paste0("  - id: ", yaml_quote(exposure_id)),
               paste0("    path: ", yaml_quote(exposure_paths[1])),
               paste0("    type: ", exposure_type))
  }

  lines <- c(lines,
             "",
             "### OUTCOME_DIR: path to outcome summary stats.",
             "# The required file format is the same as GWAS summary file.",
             paste0("OUTCOME_DIR: ", yaml_quote(outcome_dir)),
             "",
             "GWAS_REFERENCE_GENOTYPE:",
             paste0("  path: ", yaml_quote(gwas_ref)),
             paste0("  duplicated_snp_path: ", duplicated_snp_path),
             "",
             paste0("eQTL_REFERENCE_GENOTYPE: ", yaml_quote(eqtl_ref)),
             "",
             "# Parameters for colocalization analysis",
             "COLOC_SETTING:",
             paste0("  WINDOW_SIZE_BP: ", as.integer(coloc_window_size_bp)),
             paste0("  COVERAGES: ", coloc_coverages),
             paste0("  THREADS: ", as.integer(coloc_threads)),
             paste0("  CUTOFF: ", coloc_cutoff))
  writeLines(lines, config_file)
  n_qtl <- length(qtl_paths)
  n_gwas <- length(exposure_paths)
  n_outcome <- length(outcome_files)
  leo.basic::leo_log(
    "csMR step2 config written: {config_file} ({n_qtl} eQTL, {n_gwas} GWAS, {n_outcome} outcome)",
    level = "success"
  )
  normalizePath(config_file, winslash = "/", mustWork = FALSE)
}

#' Run csMR Step 3
#'
#' @description
#' Run the csMR Snakemake workflow with `R_HOME` and `R_LIBS` explicitly set
#' from the target csMR conda environment.
#'
#' @details
#' The first four arguments (`repo_dir`, `config_file`, `jobs`, `forcerun`)
#' are the most commonly adjusted by users in routine runs.
#'
#' @param repo_dir csMR repository directory.
#' @param config_file Config file path.
#' @param jobs Number of Snakemake jobs. If `NULL`, run with bare `-j` and let
#'   Snakemake use all available cores.
#' @param forcerun Optional Snakemake targets/rules to force-run. Default
#'   `NULL` means no `--forcerun` is added.
#' @param work_flow.snakefile Snakemake workflow file path. If `NULL`, use
#'   `file.path(repo_dir, "work_flow.snakefile")`.
#' @param env_name Conda env name.
#' @param conda Path to conda.
#' @param dry_run Whether to run `--dry-run`.
#' @param log_file Optional log file path.
#' @return A list with `command` (copyable shell command), `status`, `log_file`,
#'   `r_lib`, and `r_home`. Output streams to console in real time and is also
#'   saved to `log_file`.
#' @examples
#' \dontrun{
#' csMR_step3_run(repo_dir = "~/Project/software/csMR", jobs = NULL, dry_run = TRUE)
#' }
#' @export
csMR_step3_run <- function(
    repo_dir = "~/Project/software/csMR",
    config_file = "./output/csMR/step2_config/config.yml",
    jobs = NULL,
    forcerun = NULL,
    work_flow.snakefile = NULL,
    env_name = "csMR",
    conda = Sys.which("conda"),
    dry_run = FALSE,
    log_file = NULL
) {
  leo.basic::leo_log("Running csMR Step 3 ..."); t0 <- Sys.time()
  repo_dir <- normalizePath(path.expand(repo_dir), winslash = "/", mustWork = TRUE)
  if (is.null(work_flow.snakefile)) work_flow.snakefile <- file.path(repo_dir, "work_flow.snakefile")
  work_flow.snakefile <- normalizePath(path.expand(work_flow.snakefile), winslash = "/", mustWork = TRUE)
  conda <- .csmr_conda(conda)
  config_file <- normalizePath(path.expand(config_file), winslash = "/", mustWork = TRUE)
  env_prefix <- .csmr_env_prefix(conda, env_name)
  if (!nzchar(env_prefix)) stop("Cannot find conda env: ", env_name)
  r_lib <- file.path(env_prefix, "lib", "R", "library")
  r_home <- file.path(env_prefix, "lib", "R")
  env_vars <- c(paste0("R_HOME=", r_home), paste0("R_LIBS=", r_lib, ":", Sys.getenv("R_LIBS")))
  parsed_cfg <- .csmr_parse_config(config_file)
  # Validate OUTCOME_DIR: Snakemake's os.listdir() treats EVERY item as an outcome wildcard
  if (nzchar(parsed_cfg$outcome_dir) && dir.exists(parsed_cfg$outcome_dir)) {
    # Auto-remove .DS_Store (macOS Finder artifact, can reappear at any time)
    ds_store <- file.path(parsed_cfg$outcome_dir, ".DS_Store")
    if (file.exists(ds_store)) {
      unlink(ds_store)
      leo.basic::leo_log("Removed .DS_Store from OUTCOME_DIR.", level = "info")
    }
    # Hard stop if any non-.ma items remain (subdirectories, exposure files, etc.)
    all_items <- list.files(parsed_cfg$outcome_dir, all.files = FALSE)
    ma_items  <- list.files(parsed_cfg$outcome_dir, pattern = "\\.ma$")
    non_ma <- setdiff(all_items, ma_items)
    if (length(non_ma) > 0L) {
      stop(
        "OUTCOME_DIR contains non-.ma items: ", paste(non_ma, collapse = ", "),
        "\nSnakemake os.listdir() will treat ALL items as outcome wildcards. ",
        "Remove them or use a clean directory with ONLY .ma files.\n",
        "OUTCOME_DIR: ", parsed_cfg$outcome_dir
      )
    }
    if (length(ma_items) == 0L) stop("OUTCOME_DIR contains no .ma files: ", parsed_cfg$outcome_dir)
  }
  log_dir <- file.path(parsed_cfg$base_output_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(log_file)) {
    log_file <- file.path(log_dir, paste0("leo_csMR_step3_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  }
  log_file <- normalizePath(log_file, winslash = "/", mustWork = FALSE)

  # Helpful monitor hints for a separate terminal while Snakemake is running.
  tail_cmd <- paste("tail -n 80 -f", shQuote(log_file))
  leo.basic::leo_log("To monitor live logs in another terminal, run:", level = "info")
  message(tail_cmd)
  if (nzchar(Sys.which("watch"))) {
    watch_cmd <- paste("watch -n 2", shQuote(paste("tail -n 80", shQuote(log_file))))
    leo.basic::leo_log("Alternative (refresh every 2s):", level = "info")
    message(watch_cmd)
  }

  run_cmd <- function(args, log_file, append = FALSE) {
    tee_mode <- if (append) "tee -a" else "tee"
    script <- tempfile(fileext = ".sh")
    on.exit(unlink(script), add = TRUE)
    writeLines(c(
      "#!/usr/bin/env bash",
      "set -o pipefail",
      paste0("export ", env_vars),
      paste0(shQuote(conda), " ", paste(shQuote(args), collapse = " "),
             " 2>&1 | ", tee_mode, " ", shQuote(log_file))
    ), script)
    Sys.chmod(script, "0755")
    status <- system(script)
    out <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
    list(status = as.integer(status), out = out)
  }
  args <- c(
    "run", "--no-capture-output", "-n", env_name,
    "snakemake", "--directory", repo_dir,
    "-s", work_flow.snakefile,
    "--configfile", config_file
  )
  if (is.null(jobs)) {
    args <- c(args, "-j")
  } else {
    jobs <- as.integer(jobs)
    if (!is.finite(jobs) || jobs < 1L) stop("`jobs` must be a positive integer.")
    args <- c(args, "-j", as.character(jobs))
  }
  if (!is.null(forcerun)) {
    forcerun <- as.character(forcerun)
    forcerun <- unlist(strsplit(forcerun, "[,[:space:]]+"))
    forcerun <- trimws(forcerun)
    forcerun <- forcerun[nzchar(forcerun)]
    if (length(forcerun) > 0L) {
      args <- c(args, "--forcerun", forcerun)
    }
  }
  if (isTRUE(dry_run)) args <- c(args, "--dry-run")
  raw_cmd <- paste(paste(env_vars, collapse = " "), conda, paste(args, collapse = " "))
  leo.basic::leo_log("Copyable command:", level = "info")
  message(raw_cmd)
  leo.basic::leo_log("Log file: {log_file}", level = "info")
  fit <- run_cmd(args, log_file)
  if (!identical(fit$status, 0L) && any(grepl("LockException|Directory cannot be locked", fit$out, ignore.case = TRUE))) {
    leo.basic::leo_log("Detected Snakemake lock. Run `--unlock` and retry once.", level = "warning")
    unlock_fit <- run_cmd(c(
      "run", "--no-capture-output", "-n", env_name,
      "snakemake", "--directory", repo_dir,
      "-s", work_flow.snakefile,
      "--configfile", config_file,
      "--unlock"
    ), log_file, append = TRUE)
    if (!identical(unlock_fit$status, 0L)) stop("Snakemake unlock failed. Check log: ", log_file)
    fit <- run_cmd(args, log_file, append = TRUE)
  }
  if (!identical(fit$status, 0L) && any(grepl("IncompleteFilesException", fit$out, ignore.case = TRUE))) {
    leo.basic::leo_log("Detected incomplete files. Retry with `--rerun-incomplete` once.", level = "warning")
    rerun_args <- c(args, "--rerun-incomplete")
    fit <- run_cmd(rerun_args, log_file, append = TRUE)
    args <- rerun_args
  }
  if (!identical(fit$status, 0L) && any(grepl("LockException|Directory cannot be locked", fit$out, ignore.case = TRUE))) {
    leo.basic::leo_log("Detected Snakemake lock after retry. Run `--unlock` and retry once.", level = "warning")
    unlock_fit <- run_cmd(c(
      "run", "--no-capture-output", "-n", env_name,
      "snakemake", "--directory", repo_dir,
      "-s", work_flow.snakefile,
      "--configfile", config_file,
      "--unlock"
    ), log_file, append = TRUE)
    if (!identical(unlock_fit$status, 0L)) stop("Snakemake unlock failed. Check log: ", log_file)
    fit <- run_cmd(args, log_file, append = TRUE)
  }
  if (!identical(fit$status, 0L)) stop("csMR step 3 failed. Check log: ", log_file)
  leo.basic::leo_log("csMR step 3 completed successfully.", level = "success")
  leo.basic::leo_time_elapsed(t0)
  list(command = raw_cmd, status = fit$status,
       log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
       r_lib = normalizePath(r_lib, winslash = "/", mustWork = FALSE),
       r_home = normalizePath(r_home, winslash = "/", mustWork = FALSE))
}
