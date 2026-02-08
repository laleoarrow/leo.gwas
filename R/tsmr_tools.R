# F1
#' Perform LD Clumping Locally or via Reference Panel
#'
#' This function performs LD clumping using either a local PLINK binary file (\code{bfile})
#' or a 1000 Genomes super population panel. Requires the \code{ieugwasr} and \code{plinkbinr} packages.
#'
#' @param dat Data frame with columns `SNP`, `pval.exposure`; optional `id.exposure`.
#' @param pop 1000G super-pop for online clumping (AFR/AMR/EAS/EUR/SAS). Used when `bfile` is NULL.
#' @param bfile PLINK reference panel prefix for local clumping (without .bed/.bim/.fam). If set, `pop` is ignored.
#' @param plink_bin Path to PLINK executable (auto-detect via plinkbinr if NULL and `bfile` is set).
#' @param clump_kb Window size in kb. Default 10000.
#' @param clump_r2 r^2 threshold. Default 0.001.
#' @keywords tsmr:
#' @export
#' @examples
#' \dontrun{
#' # Reference:
#' # - https://github.com/MRCIEU/TwoSampleMR/issues/173
#' # - https://blog.csdn.net/xiaozheng1213/article/details/126269969
#' library(ieugwasr)
#' library(plinkbinr) # devtools::install_github("explodecomputer/plinkbinr")
#' plinkbinr::get_plink_exe()
#'
#' # Note: after using this, please check `leo_clump` column to see if they are all TRUE
#' # If it's contains F, it means no SNPs remained after clumping or something bad happened
#' }
clump_data_local <- function(dat, pop = NULL, bfile = NULL,
                             clump_kb = 10000,
                             clump_r2 = 0.001,
                             plink_bin = plinkbinr::get_plink_exe()) {
  # require(ieugwasr); require(plinkbinr)
  if (is.null(pop) & is.null(bfile)) { stop("Must indicate LD method") }
  if (!is.null(pop) & is.null(bfile)) { leo_log("Performing LD online") }
  if (is.null(pop) & !is.null(bfile)) { leo_log("Performing LD locally") }
  # if (!is.null(pop) & !is.null(bfile)) { stop("Only one LD method can be used") }

  id0 <- if ("id.exposure" %in% names(dat)) unique(dat$id.exposure)[1] else NA
  dat1 <- tryCatch(
    ieugwasr::ld_clump(dat = dplyr::tibble(rsid=dat$SNP,
                                           pval=dat$pval.exposure,
                                           id=dat$id.exposure),
                       clump_kb = clump_kb,
                       clump_r2 = clump_r2,
                       pop = pop,
                       bfile = bfile,
                       plink_bin = plink_bin),
    error = function(e) {
      if (is.na(id0) || length(id0) == 0) id0 <- "NA"
      leo_log(sprintf("ld_clump failed (id=%s): %s", id0, conditionMessage(e)), level = "danger")
      NULL
    }
  )

  if (is.null(dat1) || nrow(dat1) == 0) {
    dat$leo_clump <- FALSE
    leo_log(sprintf("No SNPs remained after clumping (id=%s)", id0), level = "danger")
    return(dat)
  }

  dat2 <- subset(dat, SNP %in% dat1$rsid)
  dat2$leo_clump <- TRUE
  return(dat2)
}

# F2
#' Extract instruments locally for MR Analysis
#'
#' Filters SNPs by p-value threshold and performs LD clumping with flexible column mapping.
#'
#' @param dat Dataframe containing GWAS summary statistics; change it to dataframe if it is data.table.
#' @param p P-value cutoff (default = 5e-8)
#' @param pop Super-population code for LD reference (default = NULL, as we recommend using loca bfile)
#' @param phenotype_col Column name for phenotype (default = "Phenotype")
#' @param snp_col Column name for SNP IDs (default = "SNP")
#' @param chr_col Column name for chromosome (default = "CHR")
#' @param pos_col Column name for position (default = "POS")
#' @param effect_allele_col Column name for effect allele (default = "A1")
#' @param other_allele_col Column name for non-effect allele (default = "A2")
#' @param beta_col Column name for effect size (default = "BETA")
#' @param se_col Column name for standard error (default = "SE")
#' @param pval_col Column name for p-values (default = "P")
#' @param eaf_col Column name for effect allele frequency (default = "EAF")
#' @param N Column name for sample size (default = "Neff", set NULL to exclude)
#' @param bfile Path to PLINK binary reference panel (default = NULL)
#' @param plink_bin Path to PLINK executable (default = NULL)
#'
#' @return Formatted exposure data ready for MR analysis
#' @export
#' @importFrom TwoSampleMR format_data
#' @examples
#' \dontrun{
#' # Custom column names example
#' extract_instruments_local(
#'   dat = gwas_data,
#'   phenotype_col = "Trait",
#'   snp_col = "rsID",
#'   chr_col = "Chromosome",
#'   pos_col = "Position"
#' )
#' }
extract_instruments_local <- function(dat, p = 5e-08, pop = NULL,
                                      phenotype_col = "Phenotype",
                                      snp_col = "SNP",
                                      chr_col = "CHR",
                                      pos_col = "POS",
                                      effect_allele_col = "A1",
                                      other_allele_col = "A2",
                                      beta_col = "BETA",
                                      se_col = "SE",
                                      pval_col = "P",
                                      eaf_col = "EAF",
                                      N = "Neff",
                                      bfile = NULL,
                                      plink_bin = NULL) {
  # check cols ---------------------------------------------------------------
  required_cols <- c(snp_col, chr_col, pos_col, effect_allele_col,
                     other_allele_col, beta_col, se_col, pval_col, eaf_col, phenotype_col)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    leo_log("Missing required columns: ", paste(missing_cols, collapse = ", "))
    stop()
  }

  if (!is.null(N) && !N %in% names(dat)) {
    stop("Sample size column '", N, "' not found in dataset")
  }

  # filter based on P ---------------------------------------------------------------
  leo_log(" - Filtering SNPs using provided P-value threshold")
  instruments <- dat[dat[[pval_col]] < p, ]
  if (nrow(instruments) == 0) {
    message("No SNPs passed the P-value threshold (p < ", p, ")")
    return(NULL)
  }
  n_iv <- nrow(instruments)
  cli::cli_alert_info(" - < { n_iv } > SNP passed the P-value threshold")

  # build args --------------------------------------------------------
  format_args <- list(
    dat = instruments,
    type = "exposure",
    phenotype_col = phenotype_col,
    snp_col = snp_col,
    chr_col = chr_col,
    pos_col = pos_col,
    effect_allele_col = effect_allele_col,
    other_allele_col = other_allele_col,
    eaf_col = eaf_col,
    beta_col = beta_col,
    se_col = se_col,
    pval_col = pval_col,
    min_pval = 1e-200
  )
  # optional
  if (!is.null(eaf_col) && eaf_col %in% names(dat)) {
    format_args$eaf_col <- eaf_col
  } else {
    message("Note: EAF column not found, harmonization may be affected")
  }

  if (!is.null(N)) {
    format_args$samplesize_col <- N
    message("Using sample size column: ", N)
  }

 instruments <- do.call(TwoSampleMR::format_data, format_args)
 # instruments <- TwoSampleMR::format_data(
 #  instruments,
 #  type = "exposure",
 #  phenotype_col = phenotype_col,
 #  snp_col = "SNP",
 #  chr_col = "CHR",
 #  pos_col = "POS",
 #  effect_allele_col = "A1",
 #  other_allele_col = "A2",
 #  eaf_col = "EAF",
 #  beta_col = "BETA",
 #  se_col = "SE",
 #  pval_col = "P",
 #  samplesize_col = N,
 #  min_pval = 1e-200
 # )
 instruments <- clump_data_local(instruments, pop = pop, bfile = bfile) # F2中执行了F1（本地Clump）
 return(instruments)
}

# F3
#' format outcome data
#'
#' @param dat a dataframe for outcome with SNP, CHR, POS, A1, A2, EAF, BETA, SE, P, Phenotype, samplesize columns
#' @param snp a str vector out of iv$SNP
#' @param N N column name for sample size (effective or observed)
#' @keywords tsmr:
#' @return a tsmr format outcome dataframe
format_outcome <- function(dat, snp = iv$SNP, N = "Neff") {
 leo_message("You should check samplesize column (effective or observed)")
 leo_message(paste0("Now using <", N, "> for samplesize indice !!!!!!!!"))
 out <- format_data(
  dat,
  snps = snp,
  type = "outcome",
  # defaut col names
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "POS",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  samplesize_col = N,
  min_pval = 1e-200
 )
 return(out)
}

# F4
#' find_proxy
#'
#' @description
#' `find_proxy` finds the proxy snp for the miss iv
#'
#' @param miss_iv iv datafram in tsmr exposure format, which can not locate snp in the outcome
#' @param miss_snp snp in the miss_iv; can be inferred using miss_iv$SNP
#' @param outcome_snp a str vector containing all snp in the outcome; this NOT the entire outcome gwas summary statistics
#' @param pop reference panel from 1kg (LDlinkR param)
#' @param gb genome build (LDlinkR param)
#' @param proxy_output_path a full path to save the proxy file when using ldlink
#' @param proxy_file pre-calculated proxy file path (full); do provide this if the proxy file is already generated !!!
#' @param token token of `LDlinkR`
#' @keywords tsmr:
#' @return a updated missiv with `proxy.snp` `proxy.effect.allele` `proxy.other.allele` `r2` col
#' @examples
#' # This function can be used when many iv can not locate corresponding snp in the outcome
#' # in tsmr analysis
#' \dontrun{
#' # iv is estracted iv via tsmr package;dat_h is a standard output of harmonise_data()
#' miss_iv <- iv[!iv$SNP %in% dat_h$SNP,]
#' miss_snp <- miss_iv$SNP
#' outcome_snp <- iri_nc$SNP
#' proxy_output_path <- "Full path to where you wanna store the LDlinkR output"
#'  proxy_iv <- find_proxy(miss_iv, miss_snp, outcome_snp,
#'                         proxy_file = "./combined_query_snp_list_grch38.txt",
#'                         proxy_output_path = NULL)
#'  # bak
#'  proxy_iv$target.snp <- proxy_iv$SNP # target snp
#'  proxy_iv$target.A1 <- proxy_iv$effect_allele.exposure
#'  proxy_iv$target.A2 <- proxy_iv$other_allele.exposure
#'  # replace for tsmr
#'  proxy_iv$SNP <- proxy_iv$proxy.snp
#'  proxy_iv$effect_allele.exposure <- proxy_iv$proxy.A1
#'  proxy_iv$other_allele.exposure <- proxy_iv$proxy.A2
#'  iv_f <- bind_rows(non_miss_iv, proxy_iv) # f for final
#'  dat_h_proxy <- harmonise_data(iv_f, out_nc_proxy)
#'  mr(dat_h_proxy) # nailed it!
#' }
#' 
#' @importFrom stats start end
#' @importFrom TwoSampleMR mr harmonise_data mr_scatter_plot 
#' @importFrom TwoSampleMR mr_forest_plot mr_funnel_plot mr_heterogeneity mr_leaveoneout mr_pleiotropy_test mr_singlesnp
find_proxy <- function(miss_iv, miss_snp, outcome_snp, proxy_file=NULL, proxy_output_path=NULL, pop="EUR", gb="grch38", token="") {
 # require(LDlinkR); require(tidyverse); require(data.table)
 # sub-part 1: using ldlink to find proxy ---------
 # this part can be skipped if the file can be provided in advance.
 if (is.null(proxy_file)) {
  # if you did not pre-calculate the proxy file, then calculate it here
  if(is.null(proxy_output_path)) {stop("Need to input where LDlink should generate the proxy file.")}
  wd <- getwd() # record work path
  setwd(proxy_output_path)
  message(paste0("Save the proxy file into path >>> ", proxy_output_path))
  message(paste0(" --- using parameter: ", pop, " (ref population)"))
  message(paste0(" --- using parameter: ", gb, " (genome_build)"))
  message(paste0("##### LDlinkR could take a while if there is many miss snps #####"))
  p <- LDlinkR::LDproxy_batch(miss_snp, # LDproxy_batch needs to write results in one txt file
                              pop = pop,  # pop = c("CEU", "TSI", "FIN", "GBR", "IBS"),
                              r2d = "r2",
                              token = token,
                              append = T, # one txt file (T) or plural txt files for each snp (F)
                              genome_build = gd)
  setwd(wd) # back to the work path
  proxy <- fread(file.path(proxy_output_path, paste0("combined_query_snp_list_", gb, ".txt"))) %>%
   filter(R2 > 0.6) %>%
   # Here, `query_A1` corresponds to `proxy_A1`; `query_A2` to `proxy_A2`
   tidyr::separate(Correlated_Alleles, sep = "[,=]", remove = FALSE, into = c("query_A1", "proxy_A1", "query_A2", "proxy_A2")) %>%
   select(query_snp, RS_Number, R2,query_A1, proxy_A1, query_A2, proxy_A2) %>%
   mutate(SNP_in_outcome = RS_Number %in% outcome_snp) %>%
   filter(SNP_in_outcome == T)
  message(paste0(" - A total of <",dim(proxy)[1], "> potential proxy snp"))
  message(paste0(" - A total of <",length(unique(proxy$query_snp)), "> miss_iv have proxy snp"))
 } else {
  proxy <- fread(proxy_file) %>% # e.g., fread("/Users/leoarrow/project/iridocyclitis/output/tsmr//combined_query_snp_list_grch38.txt")
   filter(R2 > 0.6) %>%
   tidyr::separate(Correlated_Alleles, sep = "[,=]", remove = FALSE, into = c("query_A1", "proxy_A1", "query_A2", "proxy_A2")) %>%
   select(query_snp, RS_Number, R2,query_A1, proxy_A1, query_A2, proxy_A2) %>%
   mutate(SNP_in_outcome = RS_Number %in% outcome_snp) %>%
   filter(SNP_in_outcome == T)
  message(paste0(" - A total of <",dim(proxy)[1], "> potential proxy snp"))
  message(paste0(" - A total of <",length(unique(proxy$query_snp)), "> miss_iv have proxy snp"))
 }

 # sub-part 2: match proxy ---------
 for (snp in miss_iv$SNP) { # for loop for snp in miss_iv$SNP
  message(paste(" - Locating proxy for: ", snp))
  miss_iv_line <- miss_iv[miss_iv$SNP == snp, ]
  miss_iv_A1 <- miss_iv_line$effect_allele.exposure
  miss_iv_A2 <- miss_iv_line$other_allele.exposure
  tmp_proxy <- proxy[proxy$query_snp == snp,] %>% dplyr::arrange(-R2, .by_group = T)
  if (nrow(tmp_proxy) == 0) {message(paste(" - Could not find proxy for: ", snp)); next}
  # For potential proxy snps
  match <- F
  i <- 1 # i for index
  while (!match)
  { # while loop for proxy snp
   # snp check
   proxy_snp_line <- tmp_proxy[i,]
   query_snp <- proxy_snp_line$query_snp
   proxy_snp <- proxy_snp_line$RS_Number
   if (query_snp == proxy_snp) {i <- i+1; next} # just in case
   # allele check
   query_A1 <- proxy_snp_line$query_A1
   query_A2 <- proxy_snp_line$query_A2
   proxy_A1 <- proxy_snp_line$proxy_A1
   proxy_A2 <- proxy_snp_line$proxy_A2
   logi_allele_match1 <- (miss_iv_A1 == query_A1 & miss_iv_A2 == query_A2) # condition 1
   logi_allele_match2 <- (miss_iv_A1 == query_A2 & miss_iv_A2 == query_A1) # condition 2
   if (logi_allele_match1 || logi_allele_match2) {
    miss_iv$proxy.snp[miss_iv$SNP == snp] <- proxy_snp
    miss_iv$proxy.r2[miss_iv$SNP == snp] <- proxy_snp_line$R2
    # if matched
    if (logi_allele_match1) {
     miss_iv$proxy.A1[miss_iv$SNP == snp] <- proxy_A1
     miss_iv$proxy.A2[miss_iv$SNP == snp] <- proxy_A2
    } else {
     miss_iv$proxy.A1[miss_iv$SNP == snp] <- proxy_A2
     miss_iv$proxy.A2[miss_iv$SNP == snp] <- proxy_A1
    }
    match <- T
    # if not matched
   } else {
    i <- i+1
    if (i > nrow(tmp_proxy)) {
     message(paste(" - Could not find proxy for: ", snp))
     break
    }
   }
  # proxy while loop end here
  }
  if (!match) {message(paste(" - No suitable proxy found for: ", snp))}
 # snp for loop end here
 }
 message(" ######## Done #########")
 miss_iv_with_proxy <- miss_iv %>% filter(!is.na(proxy.snp))
 located_proxy_length <- length(miss_iv_with_proxy$proxy.snp)
 miss_iv_length <- length(miss_iv$SNP)
 message(paste0(" - A total of <", located_proxy_length, "> proxy located."))
 message(paste0(" - A total of <", miss_iv_length-located_proxy_length, "> proxy not located."))
 return(miss_iv_with_proxy)
}


#'  One-Click Perform 2SMR
#'
#' `perform_mr_for_one_pair` perform a MR for one pair of exp and out
#'
#' @param dat_h harmonized data
#' @param exp a str indicating the exposure in dat_h
#' @param out a str indicating the outcome in dat_h
#' @param res_dir dir path where the result of MR analysis stored
#' @param fig_dir dir path where the figure of MR analysis stored
#' @param save_plot only save the plot if TRUE, default TRUE.
#'
#' @keywords tsmr:
#' @export
#' @return list(res_pair=res_pair, res_pair_presso=res_pair_presso)
#' @examples
#' \dontrun{
#' # dat_h is harmonized TwoSampleMR data with columns:
#' # exposure, outcome, mr_keep, beta.exposure, beta.outcome, etc.
#' out <- mr_one_pair(
#'   dat_h = dat_h,
#'   exp = "BMI",
#'   out = "T1D",
#'   save_plot = FALSE,
#'   res_dir = "./output/tsmr",
#'   fig_dir = "./figure/tsmr"
#' )
#' names(out)
#' }
#' @importFrom TwoSampleMR mr harmonise_data mr_scatter_plot mr_forest_plot mr_funnel_plot mr_heterogeneity mr_leaveoneout mr_pleiotropy_test mr_singlesnp mr_egger_regression mr_egger_regression_bootstrap default_parameters mr_leaveoneout_plot
mr_one_pair <- function(dat_h, exp = "", out = "", save_plot = TRUE, res_dir= "./output/tsmr", fig_dir="./figure/tsmr") {
  # Initialize
  pairname <- paste0(exp, " VS ", out); leo_message(paste0(" - MR Pair: ", pairname))
  if (!dir.exists(res_dir)) {dir.create(res_dir, recursive = T)}; leo_message(paste0(" - Setting Check: `res_dir` using ", res_dir))
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive = T)}; leo_message(paste0(" - Setting Check: `fig_dir` using ", fig_dir))

  # MR for one specific pair of exposure and outcome
  dat_h_pair <- dat_h %>% filter(exposure == exp, outcome == out)
  dat_h_pair <- subset(dat_h_pair, mr_keep); rownames(dat_h_pair) <- NULL
  if (nrow(dat_h_pair) < 2) {
    res_pair <- mr(dat_h_pair, method_list =  c("mr_wald_ratio"))
  } else {
    res_pair <- mr(dat_h_pair, method_list =  c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  }
  res_pair <- generate_odds_ratios(res_pair)
  res_pair <- res_pair %>%
    mutate(R2_total = scales::percent(sum(dat_h_pair$R2), accuracy = 0.01),
           Neff_max_exp = as.integer(max(dat_h_pair$samplesize.exposure[1])),
           Neff_max_out = as.integer(max(dat_h_pair$samplesize.outcome[1])))

  # mr_heterogeneity & pleiotropy test
  heterogeneity_pair <- mr_heterogeneity(dat_h_pair)
  pleiotropy_res <- mr_pleiotropy_test(dat_h_pair)
  presso_res <- tryCatch(expr = {run_mr_presso(dat_h_pair, SignifThreshold = 0.05, NbDistribution = ifelse(nrow(dat_h_pair)/0.05<1000, 1000, nrow(dat_h_pair)/0.05))},
                         error = function(e) {message("Error in MR-PRESSO: ", e$message);return(NULL)})

  # integrate results into `res_pair`
  res_pair$heterogeneity_IVW_Q <- ifelse(!is.null(heterogeneity_pair$Q[which(heterogeneity_pair$method == "Inverse variance weighted")]), heterogeneity_pair$Q[which(heterogeneity_pair$method == "Inverse variance weighted")], "NA")
  res_pair$heterogeneity_IVW_Q_df <- ifelse(!is.null(heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "Inverse variance weighted")]), heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "Inverse variance weighted")], "NA")
  res_pair$heterogeneity_IVW_Q_pval <- ifelse(!is.null(heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "Inverse variance weighted")]), heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "Inverse variance weighted")], "NA")
  res_pair$heterogeneity_Egger_Q <- ifelse(!is.null(heterogeneity_pair$Q[which(heterogeneity_pair$method == "MR Egger")]), heterogeneity_pair$Q[which(heterogeneity_pair$method == "MR Egger")], "NA")
  res_pair$heterogeneity_Egger_Q_df <- ifelse(!is.null(heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "MR Egger")]), heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "MR Egger")], "NA")
  res_pair$heterogeneity_Egger_Q_pval <- ifelse(!is.null(heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "MR Egger")]), heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "MR Egger")], "NA")
  res_pair$pleiotropy_Egger_intercept <- ifelse(!is.null(pleiotropy_res$egger_intercept), pleiotropy_res$egger_intercept, "NA")
  res_pair$pleiotropy_Egger_intercept_se <- ifelse(!is.null(pleiotropy_res$se), pleiotropy_res$se, "NA")
  res_pair$pleiotropy_Egger_intercept_pval <- ifelse(!is.null(pleiotropy_res$pval), pleiotropy_res$pval, "NA")

  # presso results
  if (!is.null(presso_res)) {
    res_pair_presso <- presso_res[[1]]$`Main MR results` %>% mutate(Exposure = exp, Outcome = out) %>% select(Exposure, Outcome, everything())
    outlier_detected <- ifelse(is.na(res_pair_presso$`Causal Estimate`[2]), "NO", "YES")
    if (outlier_detected == "YES") {
      res_pair_presso$outlier_detected <- outlier_detected
      res_pair_presso$Global_RSSobs <- presso_res[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
      res_pair_presso$Global_Pvalue <- presso_res[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
      res_pair_presso$Distortion_Coefficient <- presso_res[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
      res_pair_presso$Distortion_Pvalue <- presso_res[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
    } else {
      res_pair_presso$outlier_detected <- c(outlier_detected, "NA")
      res_pair_presso$Global_RSSobs <- c(presso_res[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs, NA)
      res_pair_presso$Global_Pvalue <- c(presso_res[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, NA)
      res_pair_presso$Distortion_Coefficient <- NA
      res_pair_presso$Distortion_Pvalue <- NA
    }
  } else {
    res_pair_presso <- NULL
  }

  # res_path <- file.path(res_dir, paste0("tsmr_", exp, "_VS_", out, ".tsv")); leo_message(paste0(" >>> Save result to: ", res_path))
  # vroom::vroom_write(res_pair, res_path)
  # res_path <- file.path(res_dir, paste0("tsmr_presso_", exp, "_VS_", out, ".tsv")); leo_message(paste0(" >>> Save result to: ", res_path))
  # vroom::vroom_write(res_pair_presso, res_path)

  if (save_plot) {
    # >>>>>>>>>>>>>  draw and save draw   >>>>>>>>>>>>>
    # p1: Scatter plot
    # p2: Forest plot
    # p3: Funnel plot
    # p4: LOO plot
    # For Forest: width = 7, height = 8 (T1D)
    # For Forest: width = 7, height = 10 (T2D)
    # For others: width = 7, height = 6
    # fig_dir <- "./figure/tsmr"
    # >>>>>>>>>>>>> Scatter Plot >>>>>>>>>>>>>>>>>>>
    # p1 <- mr_scatter_plot_modified(mr_results = res_pair, dat = dat_h_pair); p1[[1]]
    p1 <- mr_scatter_plot(mr_results = res_pair, dat = dat_h_pair)
    p1[[1]] <- p1[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "top",legend.direction = "vertical",
                     axis.title.x = element_text(size = 16, colour = "black"),  # X label
                     axis.title.y = element_text(size = 16, colour = "black"),  # Y label
                     axis.text = element_text(size = 14, colour = "black"), # x/y axis size
                     legend.title = element_text(size = 14, face = "italic"), # legend title
                     legend.text = element_text(size = 12),
                     legend.background = element_blank()
      )+ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3))
    p1[[1]][["layers"]][[1]][["aes_params"]]$colour <- "black" # error bar (|)
    p1[[1]][["layers"]][[1]][["aes_params"]]$alpha <- 0.75 # error bar (|)
    p1[[1]][["layers"]][[2]][["aes_params"]]$colour <- "black" # error bar (-)
    p1[[1]][["layers"]][[2]][["aes_params"]]$alpha <- 0.75 # error bar (-)
    p1[[1]][["layers"]][[3]][["aes_params"]]$colour <- "black" # snp point
    p1[[1]][["layers"]][[4]][["aes_params"]]$size <- 1 # geom_abline line
    p1[[1]][["layers"]][[4]][["aes_params"]]$alpha <- 0.75 # geom_abline line
    fig_path <- file.path(fig_dir, paste0("tsmr_scatter_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    ggsave(fig_path, plot = p1[[1]], width = 7, height = 6)
    # >>>>>>>>>>>>> mr_forest_plot / mr_funnel_plot >>>>>>>>>>>>>>>>>>>
    res_single_pair <- mr_singlesnp(dat_h_pair)
    p2 <- mr_forest_plot(res_single_pair)
    p3 <- mr_funnel_plot(res_single_pair)
    p2[[1]] <- p2[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "none",
                     axis.title.y = element_text(colour = "black"),
                     axis.title.x = element_text(colour = "black"),
                     axis.text = element_text(colour = "black"))
    p3[[1]] <- p3[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.title = element_text(size = 14, face = "italic"),
                     legend.text = element_text(size = 12),
                     legend.position = "top",legend.direction = "horizontal",
                     axis.title.y = element_text(size = 16,colour = "black"),
                     axis.title.x = element_text(size = 16,colour = "black"),
                     axis.text = element_text(size = 14,colour = "black"))
    fig_path <- file.path(fig_dir, paste0("tsmr_forrest_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    if (grepl("T2D*", exp)) {
      ggsave(fig_path, plot = p2[[1]], width = 7, height = 14)
    } else {
      ggsave(fig_path, plot = p2[[1]], width = 7, height = 8)}
    fig_path <- file.path(fig_dir, paste0("tsmr_funnel_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    ggsave(fig_path, plot = p3[[1]], width = 7, height = 6)
    # >>>>>>>>>>>>> LOO Plot >>>>>>>>>>>>>>>>>>>
    loo_pair = mr_leaveoneout(dat_h_pair)
    p4 <- mr_leaveoneout_plot(loo_pair)
    p4[[1]] <- p4[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "none",
                     axis.title.y = element_text(colour = "black"),
                     axis.title.x = element_text(colour = "black"),
                     axis.text = element_text(colour = "black"))
    fig_path <- file.path(fig_dir, paste0("tsmr_loo_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    if (grepl("T2D*", exp)) {
      ggsave(fig_path, plot = p4[[1]], width = 7, height = 14)
    } else {
      ggsave(fig_path, plot = p4[[1]], width = 7, height = 8)}
  }

  # reture `res-pair` and `presso` results as a list
  return(list(
    res_pair=res_pair %>% select(-id.exposure,-id.outcome) %>% select(exposure, outcome, everything()),
    res_pair_presso=res_pair_presso
    ))
}


#' Modified MR Scatter Plot
#'
#' The `mr_scatter_plot` in the `TwoSampleMR` package could be better.
#' This is a modified version of `mr_scatter_plot` in the `TwoSampleMR` package
#' When the slope of two method is really close, the scatter plot may not properly plot them!
#' Change the line type here mannually herein then.
#'
#' @param mr_results same as TwoSampleMR::mr_scatter_plot
#' @param dat same as TwoSampleMR::mr_scatter_plot
#'
#' @keywords tsmr:
#' @export
#'
#' @examples
#' \dontrun{
#' p1 <- mr_scatter_plot_modified(mr_results = res_pair, dat = dat_h_pair)
#' print(p1[[1]])
#' }
mr_scatter_plot_modified <- function(mr_results, dat) {
  # library(ggplot2);library(ggsci);library(TwoSampleMR);library(dplyr)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d) {
    d <- plyr::mutate(d)
    if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] &
                      id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if ("MR Egger" %in% mrres$method) {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome,
                                  d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    if ("MR Egger (bootstrap)" %in% mrres$method) {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome,
                                            d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }

    # debugging
    # print(mrres)

    p <- ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, y = beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                             colour = "black", alpha= 0.75, width = 0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure),
                              colour = "black", alpha= 0.75, height = 0) +
      ggplot2::geom_point(alpha= 0.75) +
      ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method, linetype = method),
                           show.legend = TRUE, size = 1.5, alpha = 0.8) +
      # ggplot2::scale_colour_manual(values = c("Inverse variance weighted" = "red",
      #                                         "MR Egger" = "blue",
      #                                         "Weighted median" = "green",
      #                                         "IVW radial" = "purple")) +
      ggplot2::scale_linetype_manual(values = c("Inverse variance weighted" = "solid", # change setting here!!!!!
                                                "MR Egger" = "solid",
                                                "Weighted median" = "solid",
                                                "IVW radial" = "dotted")) +
      ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::labs(colour = "MR Test", linetype = "MR Test", x = paste("SNP effect on", d$exposure[1]),
                    y = paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position = "top", legend.direction = "vertical") +
      ggplot2::theme(legend.position = "top",legend.direction = "vertical",
                     axis.title.x = element_text(size = 16, colour = "black"),  # X label
                     axis.title.y = element_text(size = 16, colour = "black"),  # Y label
                     axis.text = element_text(size = 14, colour = "black"), # x/y axis size
                     legend.title = element_text(size = 14, face = "italic"), # legend title
                     legend.text = element_text(size = 12),
                     legend.background = element_blank()
      )+ggplot2::guides(colour = ggplot2::guide_legend(ncol = 4), linetype = ggplot2::guide_legend(ncol = 4))
    return(p)
  })
  mrres
}

#' mrlap_one_pair
#'
#' `mrlap_one_pair` perform the MR-lap for one pair of exposure and outcome
#' [MR-lap](https://github.com/n-mounier/MRlap?tab=readme-ov-file)
#'
#' @param exposure_name exposure_name
#' @param outcome_name outcome_name
#' @param exposure_data exposure_data
#' @param outcome_data outcome_data
#' @param ld_path ldsc required file
#' @param hm3_path ldsc required file; tutorial use nomhc version list
#' @param log_path path to where you wanna store the ldsc log
#' @param log Logical, whether to save the log file (default: TRUE).
#' @param wd working directory.
#' @return A data frame containing the results of the MR-lap analysis.
#' @export
mrlap_one_pair <- function(exposure_data, outcome_data, exposure_name, outcome_name, ld_path, hm3_path, log_path, log = TRUE, wd = getwd()) {
  # Check if MRlap is available (it's in Suggests, not Imports)
  if (!requireNamespace("MRlap", quietly = TRUE)) {
    stop("Package 'MRlap' is required for this function but is not installed.\n",
         "Please install it from GitHub: remotes::install_github('n-mounier/MRlap')",
         call. = FALSE)
  }
  
  leo_message(paste0(" - MR-lap for: ", exposure_name, " VS ", outcome_name))
  if (log) {leo_message(paste0("log file would be stored at >>>"), log_path)}
  # locate iv_mrlap snps
  iv_mrlap <- extract_instruments_local(exposure_data, p = 5e-08, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
  iv_mrlap <- iv_mrlap %>%
    mutate(R2 = (2*(beta.exposure)^2*(eaf.exposure*(1-eaf.exposure)))/((se.exposure^2) * samplesize.exposure),
           F = R2*(samplesize.exposure-2)/(1-R2)) %>%
    filter(F > 10) %>%
    mutate(maf.exposure = ifelse(eaf.exposure<0.5, eaf.exposure, 1-eaf.exposure)) %>%
    filter(maf.exposure > 0.01)
  iv_mrlap$id.exposure <- iv_mrlap$exposure
  out_mrlap <- format_outcome(outcome_data, iv_mrlap$SNP)
  out_mrlap <- out_mrlap %>%
    mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>%
    subset(maf.outcome>0.01) %>%
    mutate(id.outcome = outcome)
  dat_h_mrlap <- harmonise_data(iv_mrlap, out_mrlap) %>% subset(mr_keep)

  # merge exposure and outcome with w_hm3.xxx.snplist (with/without MHC)
  hm3_snp.list <- fread(hm3_path) %>% select(SNP)
  iv_snp.list <- dat_h_mrlap %>% select(SNP)
  snp.list <- rbind(hm3_snp.list, iv_snp.list) %>% unique()
  exposure_data <- exposure_data %>% filter(SNP %in% snp.list$SNP) %>% select(-Neff)
  outcome_data <- outcome_data %>% filter(SNP %in% snp.list$SNP) %>% select(-Neff)

  # MR-lap
  wd <- getwd()
  if (log) {setwd(log_path)}
  A <- MRlap::MRlap(exposure = exposure_data,
             exposure_name = exposure_name,
             outcome = outcome_data,
             outcome_name = outcome_name,
             ld = ld_path,
             hm3 = hm3_path,
             save_logfiles = log,
             do_pruning = FALSE,
             user_SNPsToKeep = dat_h_mrlap$SNP
             )
  setwd(wd)

  # Results
  res_lap_pair <- data.frame(
    exposure = exposure_name,
    outcome = outcome_name,
    m_IVs = A[[1]]$m_IVs,
    observed_effect = A[[1]]$observed_effect,
    observed_se = A[[1]]$observed_effect_se,
    observed_p = A[[1]]$observed_effect_p,
    corrected_effect = A[[1]]$corrected_effect,
    corrected_effect_se = A[[1]]$corrected_effect_se,
    corrected_effect_p = A[[1]]$corrected_effect_p,
    test_difference = A[[1]]$test_difference,
    p_difference = A[[1]]$p_difference
  )
  return(res_lap_pair)
}

# Global variables for R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(
  "Outcome", "id.exposure", "id.outcome", "mr_keep", "beta.exposure", "beta.outcome",
  "se.outcome", "se.exposure", "a", "b", "method", "eaf.exposure", "samplesize.exposure",
  "R2", "maf.exposure", "eaf.outcome", "maf.outcome", "outcome", "SNP", "Neff",
  "proxy.snp", "effect_allele.exposure", "other_allele.exposure", "proxy.A1",
  "proxy.A2", "target.snp", "target.A1", "target.A2", "RS_Number", "Correlated_Alleles",
  "query_snp", "proxy_snp", "proxy_A1", "proxy_A2", "query_A1", "query_A2",
  "dat_h", "iri_nc", "CHR_BP", "DAF", "Exposure", "POS", "mr", "harmonise_data",
  "mr_scatter_plot", "mr_forest_plot", "mr_funnel_plot", "mr_heterogeneity",
  "mr_leaveoneout", "mr_pleiotropy_test", "mr_singlesnp",
  "Phenotype", "data", ".", "n_snp", "lo_ci", "up_ci", "or",
  "or_lci95", "or_uci95", "nsnp", "pval", "F_U", "Feature", "GI", "GRanges", "ID",
  "Importance", "MAF", "Occurence", "RAF", "RefSNP_id", "SNP_in_outcome", "biotype",
  "bp_end", "bp_start", "chr", "clean_symbol", "count", "cvlo", "cvm", "cvup", "description",
  "ensembl_gene_id", "entrez_id", "exposure", "external_gene_name", "feature",
  "feature1", "feature1_index", "feature2", "feature2_index", "first", "gd", "gene_id",
  "gene_name", "gene_symbol", "generate_odds_ratios", "genotyped_count",
  "hgnc_symbol", "idx", "imputed_count", "infer_version", "iteration", "iv", "lambda",
  "link_LD", "link_recomb", "ll", "locus", "mpfr", "n", "n_feat", "p_delong", "p_u_test", 
  "pos", "query_symbol", "r", "ranges", "rank_censor", "ratio", "row_number", "rsID", 
  "run_mr_presso", "s", "seqnames", "shap", "snpsById", "strand", "symbol", "top_n", "value"
))
