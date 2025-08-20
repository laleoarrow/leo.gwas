#' GWAS summary QC pipeline (chip + imputed)
#'
#' Perform QC on GWAS summary stats by merging imputed and genotyped data,
#' removing duplicates, checking consistency, matching to a reference panel,
#' filtering by DAF/F_U cutoffs, and saving results.
#'
#' @param summary_x2_p Path to imputed summary file (.assoc/.gz).
#' @param summary_x2_chip_p Path to genotyped (chip) summary file.
#' @param ref_p Path to reference panel (default: 1KG EAS file).
#' @param save_dir Output directory (default: \code{"~/Project/BD2025/data/qc/sex_split"}).
#' @param save_name_prefix Prefix for output file names.
#' @param DAF_cutoff Numeric. Max |DAF| allowed (default 0.2).
#' @param F_U_cutoff Numeric. Min F_U required (default 0.01).
#'
#' @details
#' Steps:
#' \itemize{
#'   \item Read and standardize imputed & chip data
#'   \item Remove duplicates (keep smallest P for genotyped, drop multiple imputed)
#'   \item Merge and label GI (imputed/genotyped)
#'   \item Match with reference panel and compute DAF
#'   \item Filter by cutoffs and save full / P<1e-6 subsets
#' }
#'
#' Logs at each step with \code{leo_log()} for dimension, duplication, NA, etc.
#'
#' @return No return value; writes QC results to \code{save_dir}.
#' @examples
#' \dontrun{
#' leo.gwas_qc("imp.assoc.gz", "chip.assoc")
#' }
leo.gwas_qc <- function(summary_x2_p, summary_x2_chip_p,
                        ref_p = "~/Project/BD2025/data/zuo/ref/1KG-EAS-EAF-chrposa2a1.txt.gz",
                        save_dir = "~/Project/BD2025/data/qc/sex_split",
                        save_name_prefix = "bd-ASA-41",
                        DAF_cutoff = 0.2, F_U_cutoff = 0.01) {
  summary_x2 <- fread(summary_x2_p) %>% leo.gwas::get_id()
  summary_x2_chip <- fread(summary_x2_chip_p) %>% leo.gwas::get_id()
  leo.basic::leo_log("QC Start: imputed = {paste(dim(summary_x2), collapse='x')}, chip = {paste(dim(summary_x2_chip), collapse='x')}")

  # remove duplicated SNP ID with same CHR:POS:A1:A2
  summary_x2_chip <- summary_x2_chip %>% # remove duplicated SNP in the chip file
    group_by(ID) %>%                     # asa has duplicated SNP while zh8 does not
    arrange(P) %>%
    dplyr::slice(1) %>%
    ungroup()

  # remove duplicated line (I mean totally duplicated here)
  summary_x2 <- summary_x2 %>% duplicated_SNP_lines(., "rm", dup_columns = c("SNP"), group_columns = c("SNP", "A1", "A2", "OR", "P"))
  summary_x2_chip <- summary_x2_chip %>% duplicated_SNP_lines(., "rm", dup_columns = c("SNP"), group_columns = c("SNP", "A1", "A2", "OR", "P"))
  leo.basic::leo_log(
    "QC for duplication: imputed = {paste(dim(summary_x2), collapse='x')}, chip = {paste(dim(summary_x2_chip), collapse='x')}"
  )

  # Add imputed or genotyped information ----
  summary_x2_qc <- bind_rows(summary_x2 %>% mutate(GI = "imputed", CHR_BP = paste0(CHR, ":", BP)),
                             summary_x2_chip %>% mutate(GI = "genotyped", CHR_BP = paste0(CHR, ":", BP))) %>% as.data.table()
  leo.basic::leo_log(
    "Combined genotyped and imputed: {paste(dim(summary_x2_qc), collapse='x')}"
  )

  # For P==0 ～ P = 0.5e-323 (minimum value in R)
  summary_x2_qc <- summary_x2_qc %>% mutate(P = case_when(
    P==0 ~ 0.5e-323,
    TRUE ~ P
  ))

  # check number of duplicated line now
  summary_x2_qc %>% duplicated_SNP_lines(., "get", dup_columns = c("CHR_BP")) %>% arrange(CHR, BP) %>% nrow() -> n
  leo.basic::leo_log("Check duplicated SNP after combining genotyped and imputed: {n}")

  ##########  step 1: 按CHR_BP来分组，如果组内genotyped和imputed都有，那么排除所有imputed
  summary_x2_qc <- summary_x2_qc %>%
    group_by(CHR_BP) %>%
    filter(!(any(GI == "genotyped") & GI == "imputed")) %>%
    ungroup()

  ########## step 2/3: For有重复的SNP，执行下面的qc（zh8不用2/3步）
  # 步骤 2：如果有多个 genotyped，保留 P 最小的
  # 步骤 3：如果有多个 imputed，全部排除
  qc_step23 <- function(summary_x2_qc){
    summary_x2_qc <- summary_x2_qc %>%
      group_by(CHR_BP) %>%
      mutate(genotyped_count = sum(GI == "genotyped"),
             imputed_count = sum(GI == "imputed")) %>%
      ungroup()
    summary_x2_qc_single <- summary_x2_qc %>% filter(genotyped_count == 1 | imputed_count == 1) # 已经合格，保留
    summary_x2_qc_genotyped_qc <- summary_x2_qc %>% filter(genotyped_count > 1) %>% # 多个genotyped的，保留P最小的
      group_by(CHR_BP) %>% arrange(P) %>% dplyr::slice(1) %>% ungroup()
    # summary_x2_qc_imputed_qc <- summary_qc %>% filter(imputed_count > 1) # 多个imputed的，这组都不要，所以不运行
    # 合并single的和多genotyped中P最小一个
    summary_qc <- rbind(summary_x2_qc_single, summary_x2_qc_genotyped_qc) %>%
      select(-genotyped_count, -imputed_count) %>%
      arrange(CHR, BP)
    return(summary_qc)
  }
  summary_qc <- qc_step23(summary_x2_qc)
  leo.basic::leo_log(
    "After 3 step QC:
   * summary_qc: dim = {paste(dim(summary_qc), collapse='x')};
   * GI table = {paste(names(table(summary_qc$GI)), as.integer(table(summary_qc$GI)), collapse=', ')};
   * GI proportion = {paste(names(prop.table(table(summary_qc$GI))), sprintf('%.2f%%', 100*as.numeric(prop.table(table(summary_qc$GI)))), collapse=', ')}"
  )

  # double check
  leo.basic::leo_log(
    "Double check:
   * duplicated SNP lines = {summary_qc %>% duplicated_SNP_lines(., 'get', dup_columns = c('CHR_BP')) %>% nrow()};
   * indel count = {summary_qc %>% fetch_indel() %>% nrow()};
   * any NA:"
  )
  summary_qc %>% any_na() %>% print()

  ######### Match MAF with reference panel (EAS 1KG) ----
  # zuo_ref <- fread("/Users/leoarrow/Project/BD2025/data/zuo/ref/1KG-EAS-EAF.txt.gz") # input with a SNP column denoting the rsid
  # zuo_ref <- zuo_ref %>% select(-CHR)
  # zuo_ref <- add_chrpos(zuo_ref, snp_col = "SNP", ref = "GRCh37") %>% drop_na(CHR)
  # zuo_ref <- get_id(zuo_ref)
  # zuo_ref %>% fwrite("/Users/leoarrow/Project/BD2025/data/zuo/ref/1KG-EAS-EAF-chrposa2a1.txt.gz")
  zuo_ref <- fread(ref_p)
  leo.basic::leo_log("Comparing with reference panel")
  summary_qc <- summary_qc %>%
    left_join(zuo_ref %>% select(ID, RAF = MAF, rsID = SNP), by = "ID") %>%  # Reference Allele Freq
    mutate(DAF = F_U - RAF) # Different Allele Freq
  leo.basic::leo_log("* any NA:")
  summary_qc %>% any_na() %>% print()
  # hist(summary_qc$DAF, breaks = 500)

  # remove non-overlapped SNPs with ref panel
  leo.basic::leo_log("Remove non-overlapped SNPs with reference panel")
  summary_qc <- summary_qc %>% drop_na(RAF)

  # QC DAF and F_U
  leo.basic::leo_log("Retain only SNPs with abs(DAF) < {DAF_cutoff} and F_U > {F_U_cutoff} reference panel")
  summary_qc <- summary_qc %>% filter(abs(DAF) < DAF_cutoff & F_U > F_U_cutoff)
  leo.basic::leo_log(
    "After matching with reference panel: dim = {paste(dim(summary_qc), collapse='x')};
   * GI table = {paste(names(table(summary_qc$GI)), as.integer(table(summary_qc$GI)), collapse=', ')};
   * GI proportion = {paste(names(prop.table(table(summary_qc$GI))), sprintf('%.2f%%', 100*as.numeric(prop.table(table(summary_qc$GI)))), collapse=', ')}"
  )

  # save
  leo.basic::leo_log("Save summary_qc to file")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  # save_name_prefix <- "bd-ASA-41"
  today_str <- format(Sys.Date(), "%Y%m%d")
  save_name <- glue::glue("{save_name_prefix}-DAF{DAF_cutoff}-FU{F_U_cutoff}-{today_str}.assoc")

  summary_qc <- summary_qc %>%
    mutate(SNP = CHR_BP) %>%
    select(-CHR_BP)
  leo_log("Save full summary_qc to {.path {file.path(save_dir, save_name)}}", level = "success")
  summary_qc %>% vroom::vroom_write(file.path(save_dir,save_name), delim = "\t")

  save_name <- glue::glue("{save_name_prefix}-DAF{DAF_cutoff}-FU{F_U_cutoff}-{today_str}-1e6.csv")
  leo_log("Save summary_qc with P < 1e-6 to {.path {file.path(save_dir, save_name)}}", level = "success")
  summary_qc %>% filter(P < 1e-6) %>% vroom::vroom_write(file.path(save_dir,save_name), delim = ",")
}
