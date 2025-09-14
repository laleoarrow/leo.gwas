#' GWAS summary QC pipeline (chip + imputed)
#'
#' Perform QC on GWAS summary stats by merging imputed and genotyped data,
#' removing duplicates, checking consistency, matching to a reference panel,
#' filtering by DAF/F_U cutoffs, and saving results.
#'
#' @param summary_x2_p Path to imputed summary file (.assoc/.gz).
#' @param summary_x2_chip_p Path to genotyped (chip) summary file.
#' @param ref_p Path to reference panel or loaded df
#' @param save_dir Output directory (default: \code{"~/Project/BD2025/data/qc/sex_split"}).
#' @param save_name_prefix Prefix for output file names.
#' @param DAF_cutoff Numeric. Max |DAF| allowed (default 0.2).
#' @param F_U_cutoff Numeric. Min F_U required (default 0.01).
#'
#' @details
#' Major steps:
#' \itemize{
#'   \item Read and standardize imputed & chip data
#'   \item Remove duplicates (keep smallest P for genotyped, drop multiple imputed)
#'   \item Merge and label GI (imputed/genotyped)
#'   \item Match with reference panel and compute DAF
#'   \item Filter by cutoffs and save full / P<1e-6 subsets and output
#' }
#'
#' Logs at each step with \code{leo_log()} for dimension, duplication, NA, etc.
#'
#' @importFrom data.table fread as.data.table
#' @importFrom dplyr group_by ungroup arrange slice mutate filter select bind_rows left_join case_when
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom vroom vroom_write
#' @importFrom glue glue
#' @importFrom leo.basic leo_log
#'
#' @return return(summary_qc); writes QC results to \code{save_dir}.
#' @examples
#' \dontrun{
#' leo.gwas_qc("imp.assoc.gz", "chip.assoc")
#' }
leo.gwas_qc <- function(summary_x2_p, summary_x2_chip_p,
                        ref_p = "/Users/leoarrow/Project/ref/1kg_maf/zuo_ref/1KG-EAS-EAF-chrposa2a1.txt.gz",
                        save_dir = "~/Project/BD2025/data/qc/sex_split",
                        save_name_prefix = "bd-ASA-41",
                        DAF_cutoff = 0.2, F_U_cutoff = 0.01) {
  # check if all file exists
  if (!file.exists(summary_x2_p)) return(leo_log("File not found: {summary_x2_p}", level = "danger"))
  if (!file.exists(summary_x2_chip_p)) return(leo_log("File not found: {summary_x2_chip_p}", level = "danger"))
  if (class(ref_p)[1]=="character"){
    if (!file.exists(ref_p)) return(leo_log("File not found: {ref_p}", level = "danger"))
  }

  # read files
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
  if (class(ref_p)[1]=="character") {
    leo.basic::leo_log("Reading reference panel from {.path {ref_p}}")
    zuo_ref <- fread(ref_p)
  } else {
    leo.basic::leo_log("Using provided reference panel object")
    zuo_ref <- ref_p
  }

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

  return(summary_qc)
}


# helper ----

#' Helper functions for `leo.gwas_qc`
#'
#' Internal utilities used in the GWAS QC pipeline.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{is_complementary(a1st, a2nd)}}{Check if two alleles form an A/T or C/G pair.}
#'   \item{\code{fetch_indel(df, type)}}{Filter indels by allele-string length.}
#'   \item{\code{fetch_non_indel(df)}}{Keep SNPs with single-base alleles.}
#'   \item{\code{duplicated_SNP_lines(df, type, dup_columns, group_columns)}}{Get/remove duplicated SNP rows.}
#'   \item{\code{slice1_SNP_lines(df, dup_columns, group_columns)}}{Within duplicated groups, keep first row.}
#'   \item{\code{fetch_same_direcrtion(df_x2, df_lg)}}{Keep same-direction effects between datasets.}
#'   \item{\code{any_na(df)}}{Count NAs per column.}
#' }
#'
#' @section Value:
#' \itemize{
#'   \item \code{is_complementary}: logical vector.
#'   \item \code{fetch_indel}, \code{fetch_non_indel}, \code{slice1_SNP_lines}, \code{duplicated_SNP_lines("rm")}: data frame.
#'   \item \code{duplicated_SNP_lines("get")}: data frame with \code{count} column.
#'   \item \code{fetch_same_direcrtion}: filtered data frame (\code{df_x2} subset).
#'   \item \code{any_na}: named numeric vector of NA counts.
#' }
#'
#' @examples
#' # Small demo dataset
#' df <- data.frame(
#'   SNP = c("rs1","rs2","rs3","rs2"),
#'   A1  = c("A","AT","C","C"),
#'   A2  = c("T","A","G","G"),
#'   OR  = c(1.2, 0.9, 1.1, 1.1)
#' )
#'
#' # Complementary alleles
#' is_complementary("A","T")  # TRUE
#' is_complementary("A","G")  # FALSE
#'
#' # Indels and non-indels
#' fetch_indel(df, "both")
#' fetch_non_indel(df)
#'
#' # Duplicates by SNP
#' duplicated_SNP_lines(df, "get", dup_columns = "SNP")
#' slice1_SNP_lines(df, dup_columns = "SNP")
#'
#' # Same-direction effects between two datasets
#' df2 <- transform(df, OR = c(1.1, 1.3, 0.8, 1.1))
#' fetch_same_direcrtion(df, df2)
#'
#' # Count NA by column
#' any_na(df)
#'
#' @name gwas_qc_helpers
#' @keywords internal
#'
#' @importFrom data.table is.data.table as.data.table
#' @importFrom dplyr filter group_by ungroup arrange mutate bind_rows slice_head
#' @importFrom dplyr across all_of
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
NULL

#' @rdname gwas_qc_helpers
is_complementary <- function(a1st, a2nd){
  return(
    (a1st == "A" & a2nd == "T") |
      (a1st == "T" & a2nd == "A") |
      (a1st == "C" & a2nd == "G") |
      (a1st == "G" & a2nd == "C")
  )
}

#' @rdname gwas_qc_helpers
fetch_indel <- function(df, type = "both"){
  if (type == "both") {
    return(df %>% filter(nchar(A1) > 1 & nchar(A2) > 1))
  }
  if (type == "A1") {
    return(df %>% filter(nchar(A1) > 1))
  }
  if (type == "A2") {
    return(df %>% filter(nchar(A2) > 1))
  }
}

#' @rdname gwas_qc_helpers
fetch_non_indel <- function(df){
  return(df %>% filter(nchar(A1) == 1 & nchar(A2) == 1))
}

#' @rdname gwas_qc_helpers
duplicated_SNP_lines <- function(df, type = "rm", dup_columns = c("SNP"), group_columns = dup_columns) {
  if (is.data.table(df)) { df <- as.data.frame(df) }
  duplicated_index <- duplicated(df[dup_columns]) | duplicated(df[dup_columns], fromLast = TRUE)
  # get duplicated SNP line
  if (type == "get") {
    df_duplicated <- df %>%
      filter(duplicated_index) %>%
      group_by(across(all_of(group_columns))) %>%
      mutate(count = n()) %>%
      ungroup()
    message("Duplicated SNP line count table")
    print(table(df_duplicated$count))
    return(df_duplicated)
  }
  if (type == "rm") {
    df_clean <- df %>%
      filter(!duplicated_index)
    return(df_clean)
  }
}

#' @rdname gwas_qc_helpers
slice1_SNP_lines <- function(df, dup_columns = c("SNP"), group_columns = dup_columns) {
  if (is.data.table(df)) { df <- as.data.frame(df) }

  # 计算重复的索引
  duplicated_index <- duplicated(df[dup_columns]) | duplicated(df[dup_columns], fromLast = TRUE)

  # 处理重复行，按group_columns分组并只保留第一行
  df_duplicated <- df %>%
    filter(duplicated_index) %>%
    group_by(across(all_of(group_columns))) %>%
    slice_head(n = 1) %>%
    ungroup()

  # 保留非重复行并合并结果
  df_new <- bind_rows(
    df %>% filter(!duplicated_index),
    df_duplicated
  ) %>% arrange(across(all_of(dup_columns)))

  return(df_new)
}

#' @rdname gwas_qc_helpers
#' @export
fetch_same_direcrtion <- function(df_x2, df_lg){
  logi <- (df_lg$OR > 1 & df_x2$OR < 1) | (df_lg$OR < 1 & df_x2$OR > 1)
  message(paste0("The number with opposite effect size：", sum(logi)))
  return(df_x2 %>% filter(!logi)) # return x2 data
}


#' @rdname gwas_qc_helpers
#' @export
any_na <- function(df){ return(df %>% map_dbl(~sum(is.na(.)))) }


