# Internal helper: run a shell command and stop on non-zero exit status.
.csmr_run <- function(cmd, args = character(), verbose = TRUE) {
  if (isTRUE(verbose)) {
    leo.basic::leo_log("Run: {cmd} {paste(shQuote(args), collapse = ' ')}")
  }
  out <- suppressWarnings(system2(cmd, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (!identical(status, 0L)) {
    tail_out <- utils::tail(out, 20L)
    stop(
      paste0(
        "Command failed: ", cmd, " ",
        paste(shQuote(args), collapse = " "),
        "\n",
        paste(tail_out, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  out
}

# Internal helper: get env prefix path from `conda env list`.
.csmr_env_prefix <- function(conda, env_name) {
  out <- .csmr_run(conda, c("env", "list"), verbose = FALSE)
  x <- trimws(out)
  x <- x[nzchar(x) & !grepl("^#|^\\*", x)]
  hit <- x[grepl(paste0("^", env_name, "\\s"), x)]
  if (length(hit) == 0L) {
    return("")
  }
  tokens <- strsplit(hit[1], "\\s+")[[1]]
  as.character(utils::tail(tokens, 1))
}

# Internal helper: best-effort removal of a possibly half-created env.
.csmr_try_remove_env <- function(conda, env_name, verbose = TRUE) {
  out <- suppressWarnings(
    system2(
      conda,
      c("env", "remove", "-n", env_name, "--yes"),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (isTRUE(verbose) && identical(status, 0L)) {
    leo.basic::leo_log("Removed partial conda env `{env_name}` before retry.", level = "warning")
  }
  invisible(NULL)
}

# Internal helper: clone csMR repository to target directory.
.csmr_clone_repo <- function(git_bin, repo_url, ref, repo_dir, verbose = TRUE) {
  dir.create(dirname(repo_dir), recursive = TRUE, showWarnings = FALSE)
  .csmr_run(
    git_bin,
    c("clone", "--depth", "1", "--branch", ref, repo_url, repo_dir),
    verbose = verbose
  )
}

# Internal helper: safely refresh repository via temp clone then atomic replace.
.csmr_replace_repo <- function(git_bin, repo_url, ref, repo_dir, verbose = TRUE) {
  tag <- format(as.integer(Sys.time()))
  tmp_repo <- paste0(repo_dir, ".tmp.", tag)
  bak_repo <- paste0(repo_dir, ".bak.", tag)
  on.exit({
    if (dir.exists(tmp_repo)) {
      unlink(tmp_repo, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  .csmr_clone_repo(
    git_bin = git_bin,
    repo_url = repo_url,
    ref = ref,
    repo_dir = tmp_repo,
    verbose = verbose
  )

  had_old <- dir.exists(repo_dir)
  if (isTRUE(had_old) && !file.rename(repo_dir, bak_repo)) {
    stop("Failed to backup existing csMR repository before overwrite.", call. = FALSE)
  }
  if (!file.rename(tmp_repo, repo_dir)) {
    if (isTRUE(had_old) && dir.exists(bak_repo)) {
      file.rename(bak_repo, repo_dir)
    }
    stop("Failed to replace csMR repository with refreshed clone.", call. = FALSE)
  }
  if (isTRUE(had_old) && dir.exists(bak_repo)) {
    unlink(bak_repo, recursive = TRUE, force = TRUE)
  }
  invisible(repo_dir)
}

# Internal helper: install csMR-required R packages inside a conda env.
.csmr_install_r_deps <- function(conda, env_name, r_lib, verbose = TRUE) {
  script <- c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "r_lib <- normalizePath(args[[1]], winslash = '/', mustWork = FALSE)",
    "dir.create(r_lib, recursive = TRUE, showWarnings = FALSE)",
    ".libPaths(c(r_lib, .libPaths()))",
    "options(ask = FALSE)",
    "options(timeout = max(1200, getOption('timeout')))",
    "repos <- c(CRAN = 'https://cloud.r-project.org')",
    "locks <- Sys.glob(file.path(r_lib, '00LOCK*'))",
    "if (length(locks)) unlink(locks, recursive = TRUE, force = TRUE)",
    "if (!requireNamespace('remotes', quietly = TRUE)) {",
    "  install.packages('remotes', repos = repos, lib = r_lib)",
    "}",
    "install_gh <- function(pkg, repo, deps = TRUE, ref = NULL) {",
    "  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))",
    "  args <- list(repo = repo, force = TRUE, upgrade = 'never', lib = r_lib, dependencies = deps)",
    "  if (!is.null(ref)) args$ref <- ref",
    "  do.call(remotes::install_github, args)",
    "  invisible(TRUE)",
    "}",
    "install_gh('phenoscanner', 'phenoscanner/phenoscanner', deps = TRUE)",
    "if (!requireNamespace('mr.raps', quietly = TRUE)) {",
    "  ok <- FALSE",
    "  for (ref in c('v0.4.1', '0.4.1')) {",
    "    fit <- try(install_gh('mr.raps', 'qingyuanzhao/mr.raps', deps = TRUE, ref = ref), silent = TRUE)",
    "    if (!inherits(fit, 'try-error')) { ok <- TRUE; break }",
    "  }",
    "  if (!ok) install_gh('mr.raps', 'qingyuanzhao/mr.raps', deps = TRUE)",
    "}",
    "install_gh('RadialMR', 'WSpiller/RadialMR', deps = TRUE)",
    "install_gh('MRMix', 'gqi/MRMix', deps = TRUE)",
    "install_gh('MRPRESSO', 'rondolab/MR-PRESSO', deps = TRUE)",
    "install_gh('ieugwasr', 'MRCIEU/ieugwasr', deps = TRUE)",
    "install_gh('TwoSampleMR', 'MRCIEU/TwoSampleMR', deps = FALSE)",
    "locks <- Sys.glob(file.path(r_lib, '00LOCK*'))",
    "if (length(locks)) unlink(locks, recursive = TRUE, force = TRUE)",
    "need <- c('getopt','coloc','ieugwasr','TwoSampleMR','phenoscanner','mr.raps','RadialMR','MRMix','MRPRESSO')",
    "miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]",
    "if (length(miss)) stop('Missing R packages: ', paste(miss, collapse = ', '))",
    "cat('csMR R dependencies installed\\n')"
  )
  tf <- tempfile(fileext = ".R")
  on.exit(unlink(tf), add = TRUE)
  writeLines(script, tf)
  .csmr_run(
    conda,
    c("run", "-n", env_name, "Rscript", tf, r_lib),
    verbose = verbose
  )
}

# Internal helper: build a compatibility env file if solver conflicts occur.
.csmr_build_compat_env <- function(env_name) {
  x <- c(
    paste0("name: ", env_name),
    "channels:",
    "  - conda-forge",
    "  - bioconda",
    "dependencies:",
    "  - python=3.10",
    "  - pybedtools",
    "  - bedtools",
    "  - snakemake=7.18.2",
    "  - r-base=4.3",
    "  - r-getopt",
    "  - r-coloc",
    "  - r-susier",
    "  - r-pkgload",
    "  - r-data.table",
    "  - r-dplyr",
    "  - r-remotes"
  )
  tf <- tempfile(fileext = ".yml")
  writeLines(x, tf)
  tf
}

# Internal helper: create conda environment with official file and fallback file.
.csmr_create_env <- function(solver_bin, conda, env_name, env_file, verbose = TRUE) {
  create_args <- c("env", "create", "--name", env_name, "--file", env_file, "--yes")
  tryCatch(
    {
      .csmr_run(solver_bin, create_args, verbose = verbose)
      invisible(TRUE)
    },
    error = function(e) {
      if (!identical(solver_bin, conda)) {
        if (isTRUE(verbose)) {
          leo.basic::leo_log("mamba solve failed; retry with conda.", level = "warning")
        }
        .csmr_try_remove_env(conda = conda, env_name = env_name, verbose = verbose)
        tryCatch(
          {
            .csmr_run(conda, create_args, verbose = verbose)
            return(invisible(TRUE))
          },
          error = function(e2) {
            if (isTRUE(verbose)) {
              leo.basic::leo_log(
                "Primary env solve failed due channel/version conflicts; retry with a minimal compatibility env file.",
                level = "warning"
              )
            }
          }
        )
      } else if (isTRUE(verbose)) {
        leo.basic::leo_log(
          "Conda solve failed; retry with a minimal compatibility env file.",
          level = "warning"
        )
      }
      .csmr_try_remove_env(conda = conda, env_name = env_name, verbose = verbose)
      compat_env <- .csmr_build_compat_env(env_name)
      on.exit(unlink(compat_env), add = TRUE)
      compat_args <- c("env", "create", "--name", env_name, "--file", compat_env, "--yes")
      .csmr_run(conda, compat_args, verbose = verbose)
      invisible(TRUE)
    }
  )
}

#' Configure Conda Environment for csMR
#'
#' @description
#' `csMR_env()` automates environment setup for the official csMR pipeline:
#' it clones/updates the csMR repository, creates (or updates) a conda
#' environment from `envs/envpy3.yml`, installs missing R dependencies used by
#' csMR scripts, and runs basic command/package checks.
#'
#' @details
#' The method follows the csMR paper workflow (Nat Commun 2024, 15:4890):
#' (1) locus-level precomputation from GWAS and single-cell eQTL summary data,
#' (2) SuSiE-based fine-mapping and colocalization, and (3) cell-type-specific
#' MR with sensitivity analyses. This function only prepares software/runtime
#' dependencies; it does not execute the full csMR analysis.
#'
#' @param repo_dir Directory where the csMR repository is stored. You can pass
#'   either the csMR repository path itself or its parent directory (for
#'   example `~/Project/software`, resolved to `~/Project/software/csMR`).
#' @param env_name Name of the conda environment.
#' @param conda Path to the conda executable. Defaults to `Sys.which("conda")`.
#' @param repo_url GitHub URL of the csMR repository.
#' @param ref Git branch/tag to checkout for csMR source code.
#' @param use_mamba Whether to prefer `mamba` for environment solving.
#' @param overwrite Whether to remove an existing environment before recreating.
#' @param update_repo Whether to check remote csMR version and overwrite local
#'   repository with a fresh clone when local is outdated/non-reproducible.
#' @param install_plink Whether to install PLINK 1.9 in the conda environment.
#' @param verbose Whether to print setup logs.
#'
#' @return A list containing `repo_dir`, `env_name`, `env_prefix`, `env_file`,
#' `r_library`, and two convenience shell snippets (`export_r_libs`,
#' `run_snakemake_example`).
#'
#' @examples
#' \dontrun{
#' # First-time setup
#' cfg <- csMR_env(repo_dir = "~/Project/software", env_name = "csMR")
#' cfg$run_snakemake_example
#'
#' # Fast health-check for an existing environment
#' csMR_env(repo_dir = "~/Project/software", env_name = "csMR", update_repo = FALSE)
#' }
#'
#' @export
csMR_env <- function(
    repo_dir = file.path(path.expand("~"), "Project", "software"),
    env_name = "csMR",
    conda = Sys.which("conda"),
    repo_url = "https://github.com/rhhao/csMR.git",
    ref = "main",
    use_mamba = TRUE,
    overwrite = FALSE,
    update_repo = TRUE,
    install_plink = TRUE,
    verbose = TRUE
) {
  repo_dir <- path.expand(repo_dir)
  if (!nzchar(conda) || !file.exists(conda)) {
    stop("`conda` was not found. Please install Miniconda/Anaconda first.", call. = FALSE)
  }
  git_bin <- Sys.which("git")
  if (!nzchar(git_bin)) {
    stop("`git` was not found in PATH.", call. = FALSE)
  }

  mamba_bin <- Sys.which("mamba")
  solver_bin <- if (isTRUE(use_mamba) && nzchar(mamba_bin)) mamba_bin else conda
  if (isTRUE(use_mamba) && !nzchar(mamba_bin) && isTRUE(verbose)) {
    leo.basic::leo_log("`mamba` not detected, fallback to conda solver.", level = "warning")
  }

  is_csmr_repo <- function(x) {
    file.exists(file.path(x, "envs", "envpy3.yml"))
  }
  if (dir.exists(repo_dir)) {
    if (!is_csmr_repo(repo_dir)) {
      nested_repo <- file.path(repo_dir, "csMR")
      if (is_csmr_repo(nested_repo) || basename(repo_dir) != "csMR") {
        repo_dir <- nested_repo
        if (isTRUE(verbose)) {
          leo.basic::leo_log("Resolve `repo_dir` to `{repo_dir}`.", level = "info")
        }
      }
    }
  } else if (basename(repo_dir) != "csMR") {
    repo_dir <- file.path(repo_dir, "csMR")
    if (isTRUE(verbose)) {
      leo.basic::leo_log("Resolve `repo_dir` to `{repo_dir}`.", level = "info")
    }
  }

  if (!dir.exists(repo_dir)) {
    .csmr_clone_repo(
      git_bin = git_bin,
      repo_url = repo_url,
      ref = ref,
      repo_dir = repo_dir,
      verbose = verbose
    )
  } else {
    if (isTRUE(update_repo)) {
      refresh_repo <- FALSE
      if (!dir.exists(file.path(repo_dir, ".git"))) {
        refresh_repo <- TRUE
        if (isTRUE(verbose)) {
          leo.basic::leo_log("Local csMR has no git metadata; overwrite with fresh clone.", level = "warning")
        }
      } else {
        refresh_repo <- tryCatch(
          {
            .csmr_run(git_bin, c("-C", repo_dir, "fetch", "origin", ref, "--depth", "1"), verbose = verbose)
            local_head <- trimws(.csmr_run(git_bin, c("-C", repo_dir, "rev-parse", "HEAD"), verbose = FALSE)[1])
            remote_head <- trimws(.csmr_run(git_bin, c("-C", repo_dir, "rev-parse", "FETCH_HEAD"), verbose = FALSE)[1])
            !identical(local_head, remote_head)
          },
          error = function(e) {
            warning("Failed to inspect local csMR version; overwrite with fresh clone: ", conditionMessage(e), call. = FALSE)
            TRUE
          }
        )
      }
      if (isTRUE(refresh_repo)) {
        if (isTRUE(verbose)) {
          leo.basic::leo_log("Detected older local csMR; overwrite with latest `{ref}`.", level = "warning")
        }
        .csmr_replace_repo(
          git_bin = git_bin,
          repo_url = repo_url,
          ref = ref,
          repo_dir = repo_dir,
          verbose = verbose
        )
      }
    }

    env_file_existing <- file.path(repo_dir, "envs", "envpy3.yml")
    if (!file.exists(env_file_existing)) {
      stop(
        "Target `repo_dir` exists but does not look like csMR repo (missing envs/envpy3.yml).",
        call. = FALSE
      )
    }
  }

  env_file <- file.path(repo_dir, "envs", "envpy3.yml")
  if (!file.exists(env_file)) {
    stop("Cannot find csMR environment file: ", env_file, call. = FALSE)
  }

  env_exists <- nzchar(.csmr_env_prefix(conda = conda, env_name = env_name))
  if (isTRUE(env_exists) && isTRUE(overwrite)) {
    .csmr_run(conda, c("env", "remove", "-n", env_name, "--yes"), verbose = verbose)
    env_exists <- FALSE
  }

  if (isTRUE(env_exists)) {
    if (isTRUE(verbose)) {
      leo.basic::leo_log(
        "Conda env `{env_name}` already exists; skip solver. Set `overwrite = TRUE` to rebuild.",
        level = "info"
      )
    }
  } else {
    .csmr_create_env(
      solver_bin = solver_bin,
      conda = conda,
      env_name = env_name,
      env_file = env_file,
      verbose = verbose
    )
  }

  env_prefix <- .csmr_env_prefix(conda = conda, env_name = env_name)
  if (!nzchar(env_prefix)) {
    stop("Failed to determine conda env prefix for: ", env_name, call. = FALSE)
  }
  r_lib <- file.path(env_prefix, "lib", "R", "library")

  if (isTRUE(install_plink)) {
    has_plink <- TRUE
    tryCatch(
      .csmr_run(conda, c("run", "-n", env_name, "plink", "--version"), verbose = FALSE),
      error = function(e) {
        has_plink <<- FALSE
      }
    )
    if (!isTRUE(has_plink)) {
      .csmr_run(
        conda,
        c("install", "-n", env_name, "-c", "bioconda", "plink=1.90b6.21", "--yes"),
        verbose = verbose
      )
    }
  }

  .csmr_install_r_deps(conda = conda, env_name = env_name, r_lib = r_lib, verbose = verbose)
  r_check_script <- c(
    "pkgs <- c('getopt','coloc','ieugwasr','TwoSampleMR','phenoscanner','mr.raps','RadialMR','MRMix','MRPRESSO')",
    "miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]",
    "if (length(miss)) stop('Missing R packages: ', paste(miss, collapse = ', '))",
    "cat('R dependencies OK\\n')"
  )
  tf_check <- tempfile(fileext = ".R")
  on.exit(unlink(tf_check), add = TRUE)
  writeLines(r_check_script, tf_check)
  .csmr_run(
    conda,
    c("run", "-n", env_name, "Rscript", tf_check),
    verbose = FALSE
  )

  .csmr_run(conda, c("run", "-n", env_name, "snakemake", "--version"), verbose = FALSE)

  if (isTRUE(install_plink)) {
    .csmr_run(conda, c("run", "-n", env_name, "plink", "--version"), verbose = FALSE)
  }

  export_r_libs <- paste0("export R_LIBS=", r_lib, ":$R_LIBS")
  run_snakemake_example <- paste(
    "conda run -n", env_name,
    "snakemake -s work_flow.snakefile --configfile config.yml -j 4"
  )

  out <- list(
    repo_dir = repo_dir,
    env_name = env_name,
    env_prefix = env_prefix,
    env_file = env_file,
    r_library = r_lib,
    export_r_libs = export_r_libs,
    run_snakemake_example = run_snakemake_example
  )
  if (isTRUE(verbose)) {
    leo.basic::leo_log("csMR environment setup completed.", level = "success")
  }
  out
}
