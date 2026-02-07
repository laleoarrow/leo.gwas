#' DuoRank PRS (Dr.PRS)
#'
#' Dr.PRS is designed to rank the importance for PRS inputs and futher optimazation.
#' It takes two none-overlap stage plink file (bfile) for machine learning modeling:
#' * lasso [to capture linear relationship]
#' * catboost [to capture non-linear relationship]
#' It also calculated PRS with plink using traditional additive model.
#'
#' @param stage1_bfile Path to the stage1 PLINK binary files. Normally it is the one generates summary data. If only 1 source of individual data is available, you can split it into 2 non-overlap datasets.
#' @param stage2_bfile Path to the stage2 PLINK binary files, i.e., the target dataset for PRS calculation.
#' @param summary_file Path to the summary statistics file or a data frame containing SNPs, alleles, and weights. Note this summary only contains refined SNP for PRS calculation.
#' @param output_dir Directory to save the output files.
#' @param snp_col Column name in the summary statistics for SNP IDs (default: "SNP.meta").
#' @param a1_col Column name in the summary statistics for effect allele (default: "A1").
#' @param weight_col Column name in the summary statistics for weights (default: "OR_Final").
#' @param weight_type Type of weights, either "OR" (default) or "Beta". If "OR", it will calculate Beta values.
#' @param plink_bin Path to the PLINK binary executable (default: `plinkbinr::get_plink_exe()`).
#' @param divide_ratio Proportion of stage1 data in iLasso (default: 0.7).
#' @param iLasso_iteration Number of iterations for iLasso model (default: 1000).
#' @param method Integer vector indicating which methods to use: 1=PLINK, 2=CatBoost, 3=iLasso. Default is all three methods.
#' @param seed Random seed (Default: 725).
#' @param dr.optimazation Logical, whether to perform greedy search to optimize PRS based on combined rank from CatBoost and iLasso (default: TRUE).
#' @param clump Place holder for future clumping function.
#' @param clump_include_hla Place holder for future clumping function.
#' @param clump_param Place holder for future clumping function.
#'
#' @return A list with:
#'   - lasso_importance, catboost_importance, combined_rank
#'   - greedy_auc_curve (data.frame: kept_snps, auc)
#'   - final_subset (character vector of SNPs)
#'   - plink_prs_file (path), logs
#' @rdname dr.prs
#' @export
dr.prs <- function(stage1_bfile="./data/zuo/bed/SNP_for_PRS/VKH-zhonghua-for_prs_snplist",
                   stage2_bfile="./data/zuo/bed/SNP_for_PRS/VKH-ASA-for_prs_snplist",
                   summary_file="./output/part1-gwas/prs/snp_for_prs.txt",
                   output_dir="./output/part1-gwas/prs/b1t2",
                   method = c(1,2,3), dr_optimazation = TRUE,
                   snp_col="SNP.meta", a1_col="A1",
                   weight_col="OR_Final", weight_type = "OR",
                   divide_ratio = 0.7,
                   iLasso_iteration = 1000, nfolds = 10,
                   plink_bin = plinkbinr::get_plink_exe(),
                   clump = T, clump_include_hla = F, clump_param = NULL, # Yet to be implemented
                   seed = 725){
  require("data.table"); require("dplyr")
  # method switch: 1=PLINK, 2=CatBoost, 3=iLasso
  method <- sort(unique(as.integer(method))); if (length(method) == 0L) method <- 1L
  stopifnot(all(method %in% c(1,2,3)))
  do_plink   <- 1L %in% method
  do_cat     <- 2L %in% method
  do_lasso   <- 3L %in% method
  do_combine <- do_cat && do_lasso
  if (do_cat) .check_catboost(context = "Method 2 (CatBoost)")
  if (do_cat || do_lasso) .check_caret(context = "method 2/3 splitting")
  leo.basic::leo_log("Methods selected: {paste(method, collapse = ', ')} (1=PLINK, 2=CatBoost, 3=iLasso)")
  # check if all file exists
  if (!file.exists(paste0(stage1_bfile, ".bed"))) return(leo.basic::leo_log("Base .bed file not found.", level = "danger"))
  if (!file.exists(paste0(stage1_bfile, ".bim"))) return(leo.basic::leo_log("Base .bim file not found.", level = "danger"))
  if (!file.exists(paste0(stage1_bfile, ".fam"))) return(leo.basic::leo_log("Base .fam file not found.", level = "danger"))
  if (!file.exists(paste0(stage2_bfile, ".bed"))) return(leo.basic::leo_log("Target .bed file not found.", level = "danger"))
  if (!file.exists(paste0(stage2_bfile, ".bim"))) return(leo.basic::leo_log("Target .bim file not found.", level = "danger"))
  if (!file.exists(paste0(stage2_bfile, ".fam"))) return(leo.basic::leo_log("Target .fam file not found.", level = "danger"))
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE); leo.basic::leo_log("Output directory created: {output_dir}")}
  if (is.character(summary_file) & !file.exists(summary_file)) return(leo.basic::leo_log("summary_file must be a file path or a data.frame.", level = "danger"))

  # clean summary stats file
  if (is.character(summary_file)) {
    summary <- data.table::fread(summary_file)
  } else {
    summary <- summary_file
  }
  summary <- summary %>%
    transmute(SNP  = .data[[snp_col]],
              A1   = .data[[a1_col]],
              BETA = if (toupper(weight_type) == "OR") base::log(.data[[weight_col]]) else .data[[weight_col]]) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    tidyr::drop_na()
  if (any(!is.finite(summary$BETA))) return(leo.basic::leo_log("Detected infinite weights after transform; re-check!", level = "danger"))
  leo.basic::leo_log("Summary statistics loaded: {nrow(summary)} SNPs; weight_type = {weight_type}")

  # check SNP overlap with SNP in bim file
  stage1_bim <- data.table::fread(paste0(stage1_bfile, ".bim"), header = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  stage2_bim <- data.table::fread(paste0(stage2_bfile, ".bim"), header = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  snp_overlap1 <- sum(summary$SNP %in% stage1_bim$SNP)
  snp_overlap2 <- sum(summary$SNP %in% stage2_bim$SNP)
  if (snp_overlap1 != nrow(summary)) leo.basic::leo_log("Only {snp_overlap1} of {nrow(summary)} SNPs found in stage1 .bim; match_rate = {round(snp_overlap1/nrow(summary), 3)}", level = "warning")
  if (snp_overlap2 != nrow(summary)) leo.basic::leo_log("Only {snp_overlap2} of {nrow(summary)} SNPs found in stage2 .bim; match_rate = {round(snp_overlap2/nrow(summary), 3)}", level = "warning")

  # ---- export SNP file once (shared by PLINK & A1 matrix extraction) ----
  plink_score_prefix <- paste0("plink_prs.", format(Sys.time(), "%Y%m%d"))
  snp_file <- file.path(output_dir, paste0(plink_score_prefix, ".txt"))
  summary %>% data.table::fwrite(snp_file, sep = "\t", row.names = FALSE, quote = FALSE)
  leo.basic::leo_log("SNP file for downstream saved to {.path {snp_file}}")

  # calculate traditional PRS with plink
  plink_prs <- NULL; plink_prs_vis <- NULL; auc_best <- NA_real_
  if (do_plink) {
    leo.basic::leo_log("Calculating traditional PRS with PLINK...")
    plink_prs(bfile = stage2_bfile,
              snp_file = snp_file,
              output_prefix = file.path(output_dir, plink_score_prefix))
    plink_prs <- fread(file.path(output_dir, paste0(plink_score_prefix, ".profile")))
    plink_prs_vis <- prs.roc_hist_density(prs_score = plink_prs$SCORESUM,
                                          phenotype = plink_prs$PHENO)
    auc_best <- as.numeric(plink_prs_vis$roc_value)
  } else {
    leo.basic::leo_log("Skip PLINK scoring (method excludes 1).", level = "warning")
  }

  # plink_get_a1_matrix
  m1 <- NULL; m2 <- NULL
  if (do_cat || do_lasso) {
    leo.basic::leo_log("Extracting genotype matrix for PRS calculation...")
    s1_a1_matrix <- plink_get_a1_matrix(bfile = stage1_bfile, summary = snp_file,
                                        output_dir = output_dir, output_name = "stage1_a1_matrix",
                                        snp_col = "SNP", a1_col = "A1")
    s2_a1_matrix <- plink_get_a1_matrix(bfile = stage2_bfile, summary = snp_file,
                                        output_dir = output_dir, output_name = "stage2_a1_matrix",
                                        snp_col = "SNP", a1_col = "A1")
    m1 <- fread(file.path(output_dir, "stage1_a1_matrix.raw"))
    m2 <- fread(file.path(output_dir, "stage2_a1_matrix.raw"))
  } else {
    leo.basic::leo_log("Skip genotype matrix extraction (no CatBoost/iLasso requested).", level = "warning")
  }

  # build CatBoost model
  train_catboost <- NULL; catboost_stage2 <- NULL; catboost_rank <- NULL
  if (do_cat) {
    train_catboost <- catboost_prs(m1, divide = TRUE, divide_ratio = divide_ratio)
    catboost_stage2 <- catboost_prs_target(m2, model = train_catboost$model)
    leo.basic::leo_log("CatBoost PRS AUC on target set: {round(catboost_stage2$perf$roc_value, 4)}")
    catboost_rank <- catboost_prs_rank(model = train_catboost$model,
                                       pool = train_catboost$train_pool,
                                       pool_df = train_catboost$x_train)
  } else {
    leo.basic::leo_log("CatBoost disabled (method excludes 2).", level = "warning")
  }

  # build iLasso model
  train_lasso <- NULL; lasso_stage2 <- NULL; lasso_rank <- NULL
  if (do_lasso) {
    leo.basic::leo_log("Building iLasso model for PRS ranking...")
    train_lasso <- lasso_prs(m1, divide_ratio = divide_ratio, iterative = iLasso_iteration, nfolds = nfolds)
    if (is.null(train_lasso$model) || is.null(train_lasso$lambda)) {
      leo.basic::leo_log("iLasso produced no valid model; skip Lasso target scoring.", level = "warning")
      lasso_stage2 <- list(pred_df = NULL, perf = list(roc_value = NA_real_))
      lasso_rank   <- tibble::tibble(Feature = character(0), count = integer(0))
    } else {
      lasso_stage2 <- lasso_prs_target(m2, model = train_lasso$model, lambda = train_lasso$lambda, score_type = "link")
      leo.basic::leo_log("Lasso PRS AUC on target set: {round(lasso_stage2$perf$roc_value, 4)}")
      lasso_rank <- lasso_prs_rank(model = train_lasso$model, rank = train_lasso$feature_rank, auc_history = train_lasso$auc_history)
    }
    # lasso_vis1 <- gglasso(no_cv_lasso = train_lasso$model$glmnet.fit, cv_lasso = train_lasso$model)
    # lasso_vis2 <- ggcvlasso(train_lasso$model)
  } else {
    leo.basic::leo_log("iLasso disabled (method excludes 3).", level = "warning")
  }

  # ---- combine ranks if both available ----
  if (dr_optimazation) {
    rank_combined <- NULL
    id_deleted <- character(0)
    if (do_cat && do_lasso && do_plink) {
      # overall rank
      rank1 <- catboost_rank$importance_df %>% mutate(rank1 = dplyr::dense_rank(Importance))
      rank2 <- train_lasso$feature_rank %>% mutate(rank2 = dplyr::dense_rank(count))
      rank_combined <- left_join(rank1 %>% select(Feature, rank1),
                                 rank2 %>% select(Feature, rank2),
                                 by = "Feature") %>%
        mutate(rank_combined = combine_rank(rank1, rank2,
                                            auc1 = as.numeric(train_catboost$perf$roc_value),
                                            auc2 = as.numeric(train_lasso$perf_test$roc_value)),
               rank_censor = dplyr::dense_rank(rank_combined)) %>%
        arrange(rank_censor)
      # Optimizing PRS
      id_deleted <- character(0)               # chr:pos that have been deleted
      improve_ids    <- character(0)           # all chr:pos that led to improvement
      improve_labels <- character(0)

      base_dir <- file.path(output_dir, "dr.prs"); dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

      leo.basic::leo_log("---- Start Greedy search to optimize PRS based on censor rank ----")
      leo.basic::leo_log("Initial AUC = {round(auc_best,4)}; {nrow(rank_combined)} SNP{?s} to try")

      round_idx <- 0L
      repeat {
        round_idx <- round_idx + 1L
        round_dir <- file.path(base_dir, sprintf("round%d", round_idx))
        dir.create(round_dir, showWarnings = FALSE, recursive = TRUE)
        leo.basic::leo_log("Round {sprintf('%d', round_idx)} -> temp files at {.path {round_dir}}")

        improved <- FALSE
        for (j in seq_len(nrow(rank_combined))) {
          # current SNP (label vs id)
          snp_to_drop    <- rank_combined$Feature[j]                              # e.g., "22:30253222_C(/A)"
          snp_id_to_drop <- stringr::str_extract(snp_to_drop, "^[0-9]+:[0-9]+")   # e.g., "22:30253222"
          if (snp_id_to_drop %in% id_deleted) next

          rank_tag     <- rank_combined$rank_censor[j]
          rank_tag_str <- sprintf("%d.", as.integer(rank_tag))

          leo.basic::leo_log("Try delete --> {.emph {snp_to_drop}}")

          # candidate summary = original - deleted - current
          summary_i <- summary %>% dplyr::filter(!SNP %in% c(id_deleted, snp_id_to_drop))

          # file stubs under round dir
          file_stub  <- sprintf("%s%s@%s", rank_tag_str, plink_score_prefix, snp_id_to_drop)
          out_prefix <- file.path(round_dir, file_stub)
          snp_file   <- file.path(round_dir, paste0(file_stub, ".txt"))

          # write SNP file
          data.table::fwrite(summary_i, snp_file, sep = "\t", quote = FALSE, row.names = FALSE)
          leo.basic::leo_log("SNP file saved: {.path {snp_file}}")

          # run PRS & evaluate
          plink_prs(bfile = stage2_bfile, snp_file = snp_file, output_prefix = out_prefix)
          prof <- data.table::fread(paste0(out_prefix, ".profile"))
          vis  <- prs.roc_hist_density(prs_score = prof$SCORESUM, phenotype = prof$PHENO)
          auc_i <- as.numeric(vis$roc_value)

          if (auc_i >= auc_best) {
            # accept and stop this round
            auc_best      <- auc_i
            id_deleted    <- c(id_deleted, snp_id_to_drop)
            improve_ids   <- c(improve_ids, snp_id_to_drop)
            improve_labels<- c(improve_labels, snp_to_drop)
            improved      <- TRUE
            leo.basic::leo_log("ACCEPT: AUC={round(auc_i,4)}; move to next round", level = "success")
            break
          } else {
            leo.basic::leo_log("REJECT: AUC={round(auc_i,4)}", level = "danger")
          }
        }
        if (!improved) { leo.basic::leo_log("No improvement this round. Stop.", level = "warning"); break }
        leo.basic::leo_log("Round {round_idx} finished with improvement; continue.", level = "success")
      }

      # final report: who contributed improvements
      if (length(improve_ids) == 0) {
        leo.basic::leo_log("No improvement over baseline. Best AUC = {round(auc_best,4)}", level = "warning")
      } else {
        leo.basic::leo_log("Best AUC = {round(auc_best,4)}; improvements came from dropping:", level = "success")
        # concise print
        print(tibble::tibble(order = seq_along(improve_ids), snp_label = improve_labels, snp_id = improve_ids))
      }
    }
  } else {
    leo.basic::leo_log("Skip DuoRank optimization (dr_optimazation = FALSE).", level = "warning")
  }


  # return
  return(list(combined_rank = if (exists("rank_combined")) rank_combined else NULL,
              # CatBoost
              catboost = if (do_cat) list(model = train_catboost$model,
                                          stage2_pred = catboost_stage2$pred_df,
                                          stage2_perf = catboost_stage2$perf,
                                          rank  = catboost_rank) else NULL,
              # iLasso
              ilasso = if (do_lasso) list(model = train_lasso$model,
                                          lambda = train_lasso$lambda,
                                          stage2_pred = lasso_stage2$pred_df,
                                          stage2_perf = lasso_stage2$perf,
                                          rank  = lasso_rank) else NULL,
              # PLINK
              plink = if (do_plink) list(prs_file = file.path(output_dir, paste0(plink_score_prefix, ".profile")),
                                         snp_file = snp_file,
                                         prs_df   = plink_prs,
                                         perf     = plink_prs_vis) else NULL,
              final_snp = summary %>% dplyr::filter(!SNP %in% (if (exists("id_deleted")) id_deleted else character(0))) %>% dplyr::pull(SNP),
              final_best_auc = if (do_cat && do_lasso && do_plink && isTRUE(dr_optimazation) && exists("auc_best")) auc_best
                               else if (do_plink) as.numeric(plink_prs_vis$roc_value) else NA_real_,
              auc = list(catboost   = if (do_cat)   as.numeric(catboost_stage2$perf$roc_value) else NA_real_,
                         ilasso     = if (do_lasso) as.numeric(lasso_stage2$perf$roc_value)    else NA_real_,
                         full_plink = if (do_plink) as.numeric(plink_prs_vis$roc_value)        else NA_real_)
  ))
}


#' combine two rank
#' @export
#' @rdname dr.prs
combine_rank <- function(rank1, rank2, auc1 = NULL, auc2 = NULL){
  r1 <- as.numeric(rank1); r2 <- as.numeric(rank2)
  if (is.null(auc1) || is.null(auc2) || !is.finite(auc1) || !is.finite(auc2)) {
    res <- (r1 + r2) / 2
    leo.basic::leo_log("Combine mode: consensus (mean rank)", level = "success")
  } else {
    w1 <- auc1 / (auc1 + auc2); w2 <- 1 - w1
    res <- w1 * r1 + w2 * r2
    leo.basic::leo_log("Combine mode: AUC-weighted ranks; w1={round(w1,3)}, w2={round(w2,3)}", level = "success")
  }
  return(res)
}

# catboost ----
#' CatBoost PRS utilities
#'
#' Train, apply, and rank CatBoost polygenic risk score (PRS) models.
#'
#' These functions require installed \pkg{catboost}. Training with \code{divide = TRUE}
#' also requires \pkg{caret} for stratified data splitting.
#'
#' @param a1_matrix A1 matrix (with PHENOTYPE and ID columns).
#' @param divide Logical; split into train/val sets. Default FALSE.
#' @param divide_ratio Numeric; train proportion if divide=TRUE. Default 0.5.
#' @param model Trained CatBoost model.
#' @param pool CatBoost pool.
#' @param pool_df Data used to build \code{pool}.
#' @param types Character; any of "FeatureImportance","ShapValues","Interaction".
#' @param top_k Integer; top features to keep. Default NULL.
#'
#' @return Each function returns a list; contents depend on task (training, target prediction, ranking).
#'
#' @seealso \code{\link[catboost]{catboost.train}}, \code{\link[catboost]{catboost.get_feature_importance}}
#' @export
#' @rdname catboost_prs
#' @importFrom dplyr mutate select across
catboost_prs <- function(a1_matrix, divide = F, divide_ratio = 0.5){
  .check_catboost()
  .check_caret()
  df <- a1_matrix %>% mutate(PHENOTYPE = ifelse(PHENOTYPE == 1, 0L, 1L)) %>% select(-c(FID, IID, PAT, MAT, SEX))  # 0 is ctrl/ 1 is control
  if (divide) {
    set.seed(725)
    train_idx <- caret::createDataPartition(df$PHENOTYPE, p = divide_ratio, list = FALSE)
    train_df <- df[train_idx, ]
    val_df  <- df[-train_idx, ]

    x_train <- train_df %>% select(-PHENOTYPE) %>% mutate(across(everything(), as.numeric))
    y_train <- train_df$PHENOTYPE
    x_val  <- val_df  %>% select(-PHENOTYPE) %>% mutate(across(everything(), as.numeric))
    y_val  <- val_df$PHENOTYPE

    train_pool <- catboost::catboost.load_pool(data = x_train, label = y_train)
    val_pool  <- catboost::catboost.load_pool(data = x_val,  label = y_val)
  } else {
    train_df <- df
    x_train <- train_df %>% select(-PHENOTYPE) %>% mutate(across(everything(), as.numeric))
    y_train <- train_df$PHENOTYPE
    train_pool <- catboost::catboost.load_pool(data = x_train, label = y_train)
  }

  model <- catboost::catboost.train(train_pool, NULL, params = list(loss_function = "Logloss",
                                                          eval_metric   = "AUC",
                                                          iterations    = 1000,
                                                          depth         = 7,
                                                          l2_leaf_reg   = 5,
                                                          learning_rate = 0.01,
                                                          random_seed   = 725,
                                                          metric_period = 50,
                                                          # generazation
                                                          bootstrap_type = "Bernoulli",
                                                          subsample      = 0.8,
                                                          rsm            = 0.8,
                                                          train_dir      = file.path(tempdir(), "catboost_model") ) )
  if (divide) {
    pred_score <- catboost::catboost.predict(model, val_pool)
    perf <- prs.roc_hist_density(pred_score, y_val, ctrl_case_level = c(0,1))
  } else {
    pred_score <- catboost::catboost.predict(model, train_pool)
    leo.basic::leo_log("No validation set; performance evaluated on training set.", level = "warning")
    perf <- prs.roc_hist_density(pred_score, y_train, ctrl_case_level = c(0,1))
  }
  return(list(model = model, perf = perf,
              ratio = if (divide) divide_ratio else NULL,
              train_pool = train_pool,
              x_train = x_train,
              y_train = y_train,
              val_pool = if (divide) val_pool else NULL,
              x_val = if (divide) x_val else NULL,
              y_val = if (divide) y_val else NULL) )
}

#' @rdname catboost_prs
#' @export
catboost_prs_target <- function(a1_matrix, model){
  .check_catboost()
  df <- a1_matrix %>% mutate(PHENOTYPE = ifelse(PHENOTYPE == 1, 0L, 1L)) %>% select(-c(FID, IID, PAT, MAT, SEX))
  test_df <- df
  x_test <- test_df %>% select(-PHENOTYPE) %>% mutate(across(everything(), as.numeric))
  y_test <- test_df$PHENOTYPE
  test_pool  <- catboost::catboost.load_pool(data = x_test,  label = y_test)
  pred_score <- catboost::catboost.predict(model, test_pool)
  perf <- prs.roc_hist_density(pred_score, y_test, ctrl_case_level = c(0,1))

  pred_df <- data.frame(FID = a1_matrix$FID,
                        PHENOTYPE = y_test,
                        PRED_SCORE = pred_score)
  return(list(pred_df = pred_df, perf = perf,
              target_pool = test_pool,
              target_x = x_test,
              target_y = y_test) )
}

#' @rdname catboost_prs
#' @export
catboost_prs_rank <- function(model, pool, pool_df,
                              types = c("FeatureImportance", "ShapValues", "Interaction"),
                              top_k = NULL) {
  .check_catboost()
  return_list <- list(); types <- as.vector(types)
  for (type in types) {
    leo.basic::leo_log("Calculating CatBoost importance (type = {type})...")
    importance <- catboost::catboost.get_feature_importance(model,
                                                  pool = pool,
                                                  type = type)
    # FeatureImportance
    if (type == "FeatureImportance") {
      importance_df <- data.frame(Feature = rownames(importance),
                                  Importance = importance) %>% arrange(Importance)
      if (!is.null(top_k)) importance_df <- importance_df %>% top_n(top_k, wt = Importance)
      importance_plot <- importance_df %>%
        mutate(Feature = factor(Feature, levels = Feature)) %>%
        ggplot(aes(x = Importance, y = Feature)) +
        geom_bar(stat = "identity", fill = "#D5563F") +
        labs(title = "CatBoost Feature Importance", x = NULL, y = NULL) +
        theme_bw() +
        # scale_x_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
        theme(plot.title.position = "panel", # or "plot"
              plot.title = element_text(hjust = 0, size = 20, face = "bold"),
              axis.title = element_text(size = 16, color = "black"),
              axis.text = element_text(size = 13, color = "black"),
              axis.line = element_line(linewidth = .5, color = "black"),
              legend.position = c(0.85, 0.7), # Place legend inside the plot at (0.8, 0.8)
              legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
              legend.key = element_blank(), # Remove border around legend keys
              legend.text = element_text(size = 12),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank() )
      return_list$importance_df <- importance_df
      return_list$importance_rank <- importance_plot
    }
    # ShapValues
    if  (type == "ShapValues") {
      shap_mat  <- importance[, -ncol(importance)] # last column is bias/expected_value
      colnames(shap_mat) <- colnames(pool_df)

      shap_df <- as_tibble(shap_mat) %>%
        mutate(id = row_number()) %>%
        tidyr::pivot_longer(-id, names_to = "feature", values_to = "shap")
      geno_df <- as_tibble(pool_df) %>%
        mutate(id = row_number()) %>%
        tidyr::pivot_longer(-id, names_to = "feature", values_to = "geno")

      mean_abs <- colMeans(abs(shap_mat))
      top_feats <- if (is.null(top_k) || top_k >= length(mean_abs)) {
        names(mean_abs)
      } else {
        names(sort(mean_abs, decreasing = T))[1:top_k]
      }
      # plot
      ## merge & order
      plot_df <- shap_df %>% dplyr::left_join(geno_df, by = c("id","feature"))
      feat_order <- tibble::tibble(feature = names(mean_abs),
                                   mean_shap = as.numeric(mean_abs)) %>%
        dplyr::arrange(mean_shap) %>%
        dplyr::pull(feature)
      ## 1) mean SHAP
      mean_df <- tibble::tibble(feature = names(mean_abs),
                                mean_shap = as.numeric(mean_abs)) %>%
        dplyr::mutate(feature = factor(feature, levels = feat_order))
      mean_shap <- mean_df %>%
        ggplot(aes(x = mean_shap, y = feature)) +
        geom_col(fill = "steelblue") +
        labs(x = NULL, y = NULL, title = "CatBoost Mean SHAP") +
        theme_bw() +
        theme(plot.title.position = "panel",
              plot.title = element_text(hjust = 0, size = 20, face = "bold"),
              axis.title = element_text(size = 16, color = "black"),
              axis.text = element_text(size = 13, color = "black"),
              axis.line = element_line(linewidth = .5, color = "black"),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank())
      ## 2) jitter beeswarm
      jitter_shap <- plot_df %>%
        dplyr::mutate(geno = factor(geno),
                      feature = factor(feature, levels = feat_order)) %>%
        ggplot(aes(x = feature, y = shap, color = geno)) +
        geom_jitter(width = 0.25, alpha = 0.5, size = 0.6) +
        ggsci::scale_color_npg() +
        coord_flip() +
        labs(x = NULL, y = NULL, title = "CatBoost SHAP", color = "Genotype") +
        theme_bw() +
        theme(plot.title.position = "panel",
              plot.title = element_text(hjust = 0, size = 20, face = "bold"),
              axis.title = element_text(size = 16, color = "black"),
              axis.text = element_text(size = 13, color = "black"),
              axis.line = element_line(linewidth = .5, color = "black"),
              legend.position = c(0.85, 0.7),
              legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
              legend.background = element_rect(fill = "gray90"),
              legend.key = element_blank(),
              legend.text = element_text(size = 12)) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
      ## 3) quasi-beeswarm
      bee_shap <- plot_df %>%
        dplyr::mutate(geno = factor(geno),
                      feature = factor(feature, levels = feat_order)) %>%
        ggplot(aes(x = shap, y = feature, color = geno)) +
        ggbeeswarm::geom_quasirandom(width = 0.2, varwidth = T, alpha = 0.4, size = 0.7) +
        geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey60") +
        ggsci::scale_color_npg() +
        labs(x = NULL, y = NULL, title = "CatBoost SHAP Beeswarm", color = "Genotype") +
        theme_classic() +
        theme(plot.title.position = "panel",
              plot.title = element_text(hjust = 0, size = 20, face = "bold"),
              axis.title = element_text(size = 16, color = "black"),
              axis.text = element_text(size = 13, color = "black"),
              axis.line = element_line(linewidth = .5, color = "black"),
              legend.position = c(0.85, 0.7),
              legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
              legend.key = element_blank(),
              legend.text = element_text(size = 12)) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
      ## 4) genotype ratio (same y order)
      geno_ratio <- plot_df %>%
        dplyr::mutate(geno = factor(geno),
                      feature = factor(feature, levels = feat_order)) %>%
        dplyr::count(feature, geno, name = "n") %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(ratio = n / sum(n)) %>%
        dplyr::ungroup()
      ratio_bar <- geno_ratio %>%
        ggplot(aes(y = feature, x = ratio, fill = geno)) +
        geom_col(width = 0.7) +
        scale_y_discrete(limits = feat_order) +
        scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
        ggsci::scale_fill_npg() +
        labs(x = "Genotype ratio", y = NULL, fill = "Genotype", title = "Genotype Composition") +
        theme_classic() +
        theme(plot.title.position = "panel",
              plot.title = element_text(hjust = 0, size = 16, face = "bold"),
              axis.title = element_text(size = 14, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              axis.line = element_line(linewidth = .5, color = "black"),
              legend.position = "none")

      return_list$shap_matrix_df <- shap_mat
      return_list$shap_mean_df <- mean_df
      return_list$shap_mean <- mean_shap
      return_list$shap_jitter <- jitter_shap
      return_list$shap_bee <- bee_shap
      return_list$shap_ratio <- ratio_bar
    }
    # CatBoost Interaction
    if (type == "Interaction") {
      feats <- colnames(pool_df %>% dplyr::select(-dplyr::any_of("PHENOTYPE")))
      inter_df <- tibble::as_tibble(importance, .name_repair = "minimal")
      if (any(inter_df$feature1_index == 0L | inter_df$feature2_index == 0L)) {
        inter_df <- inter_df %>%
          dplyr::mutate(feature1_index = feature1_index + 1L,
                        feature2_index = feature2_index + 1L)
      }
      inter_df <- inter_df %>%
        dplyr::mutate(feature1 = feats[feature1_index],
                      feature2 = feats[feature2_index])
      # df process
      pairs <- inter_df %>%
        dplyr::mutate(a = pmin(feature1, feature2),
                      b = pmax(feature1, feature2)) %>%
        dplyr::group_by(a, b) %>%
        dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop")
      inter_sym <- dplyr::bind_rows(dplyr::transmute(pairs, feature1 = a, feature2 = b, score),
                                    dplyr::transmute(pairs, feature1 = b, feature2 = a, score)) %>%
        dplyr::mutate(feature1 = factor(feature1, levels = feats),
                      feature2 = factor(feature2, levels = feats))
      all_pairs <- tidyr::expand_grid(feature1 = feats, feature2 = feats) %>%
        dplyr::left_join(inter_sym, by = c("feature1","feature2")) %>%
        dplyr::mutate(feature1 = factor(feature1, levels = feats),
                      feature2 = factor(feature2, levels = feats))
      p_inter <- ggplot(all_pairs, aes(x = feature1, y = feature2, fill = score)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "red", na.value = "grey") +  # <- NA 用灰色
        theme_minimal() +
        labs(x = NULL, y = NULL, title = "CatBoost Feature Interaction",fill = "Score") +
        theme(
          axis.text = element_text(size = 12, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0, size = 16, face = "bold"),
          axis.title = element_text(size = 14, color = "black"),
          legend.position = "right"
        )
      return_list$interaction_df <- inter_df
      return_list$interaction_plot <- p_inter
    }
  }
  leo.basic::leo_log("Done!", level = "success")
  return(return_list)
}

# lasso ----
#' Iterative Lasso PRS (iLasso) utilities
#'
#' Train, apply and visualize iterative Lasso-based PRS models with success gating.
#'
#' @param a1_matrix data.frame; PLINK A1 matrix.
#' @param divide_ratio numeric; train fraction per iteration (default 0.7).
#' @param iterative integer; number of iterations (default 100) - If set to 1, then regular lasso is performed.
#' @param nfolds integer; CV folds for glmnet (default 10).
#' @param seed integer; base random seed (default 725).
#' @param score_type character; "link" (linear score, default) or "response" (probability).
#' @param model Trained glmnet model (for target/rank functions).
#' @param lambda Best lambda from CV.
#' @param rank Feature frequency table.
#' @param auc_history Tibble of iteration results.
#'
#' @return Each function returns a list:
#' \itemize{
#'   \item \code{lasso_prs}: model, lambda, perf_train, perf_test, success_n, attempts, auc_history, feature_rank
#'   \item \code{lasso_prs_target}: pred_df, perf, target_x, target_y
#'   \item \code{lasso_prs_rank}: plots for feature frequency and diagnostics
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link[pROC]{roc}}
#' @export
#' @rdname lasso_prs
#' @importFrom dplyr mutate select everything across bind_rows arrange desc tibble slice
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef predict wilcox.test
#' @importFrom pROC roc roc.test
lasso_prs <- function(a1_matrix, divide_ratio = 0.7, iterative = 100, nfolds = 10, score_type = "link", seed = 725){
  .check_caret()
  # prepare base frame
  base_df <- a1_matrix %>% dplyr::mutate(PHENOTYPE = ifelse(PHENOTYPE == 1, 0L, 1L)) %>% dplyr::select(-c(FID, IID, PAT, MAT, SEX))
  n_total <- nrow(base_df); p <- ncol(base_df) - 1
  leo.basic::leo_log("iLasso start: iterative={iterative}, nfolds={nfolds}, divide_ratio={divide_ratio}, n={n_total}, nfeature={p}")

  # trackers
  attempts <- 0L; success_n <- 0L
  best_model <- list(test_auc = -Inf, model = NULL, lambda = NULL, train_idx = NULL, test_idx = NULL)
  feature_count <- setNames(integer(p), colnames(base_df %>% dplyr::select(-PHENOTYPE)))
  auc_history <- dplyr::tibble(iteration = integer(0), train_auc = numeric(0), test_auc = numeric(0),
                               pass = logical(0), p_u_test = numeric(0), p_delong = numeric(0), n_feat = integer(0))

  # iterative loop --- i.e., iLasso (upgrade from previous version [Lu et al, 2024 BSPC, DOI: 10.1016/j.bspc.2024.106271])
  pb <- cli::cli_progress_along(seq_len(iterative), name = "iLasso", clear = TRUE)
  for (i in pb) {
    attempts <- attempts + 1L; set.seed(seed + i)

    # per-iteration split
    idx <- caret::createDataPartition(base_df$PHENOTYPE, p = divide_ratio, list = FALSE)
    train_df <- base_df[idx, ]; test_df <- base_df[-idx, ]

    x_train <- train_df %>% dplyr::select(-PHENOTYPE) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>% as.matrix()
    y_train <- train_df$PHENOTYPE
    x_test <- test_df %>% dplyr::select(-PHENOTYPE) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>% as.matrix()
    y_test <- test_df$PHENOTYPE

    # randomized CV folds
    cv_fit <- glmnet::cv.glmnet(x_train, y_train,
                                family = "binomial",
                                alpha = 1,
                                type.measure = "deviance",
                                nfolds = nfolds,
                                intercept = TRUE)
    model <- cv_fit$glmnet.fit
    lambda <- cv_fit$lambda.min

    # selected features
    coef_mat <- as.matrix(stats::coef(model, s = lambda))
    nonzero_idx <- which(coef_mat[-1, 1] != 0)
    if (!length(nonzero_idx)) {
      auc_history <- dplyr::bind_rows(auc_history, dplyr::tibble(iteration = i, train_auc = NA_real_, test_auc = NA_real_,
                                                                 pass = FALSE, p_u_test = NA_real_, p_delong = NA_real_, n_feat = 0L))
      if (iterative > 1L) next # keep going when iterative == 1L (we still want a valid single-run model)
    }

    # predictions
    score_train <- drop(stats::predict(cv_fit, newx = x_train, s = lambda, type = score_type))
    score_test  <- drop(stats::predict(cv_fit, newx = x_test,  s = lambda, type = score_type))
    roc_train  <- suppressMessages(pROC::roc(y_train, score_train, quiet = TRUE))
    roc_test   <- suppressMessages(pROC::roc(y_test,  score_test,  quiet = TRUE))

    # success rules: (1) test separation; (2) no train-test AUC diff
    p_u  <- stats::wilcox.test(score_test[y_test == 0], score_test[y_test == 1], exact = FALSE)$p.value # u test
    p_dl <- suppressMessages(pROC::roc.test(roc_train, roc_test, method = "delong")$p.value) # delong test
    ok <- (p_u < 0.05) && (p_dl >= 0.05)
    auc_history <- dplyr::bind_rows(auc_history, dplyr::tibble(iteration = i,
                                                               train_auc = as.numeric(roc_train$auc),
                                                               test_auc  = as.numeric(roc_test$auc),
                                                               pass = ok, p_u_test = p_u, p_delong = p_dl,
                                                               n_feat = length(nonzero_idx)))

    # single-run: always keep this model (overwrite to ensure exactly-one model)
    if (iterative == 1) {
      success_n = as.integer(ok)
      if (length(nonzero_idx)) feature_count[nonzero_idx] <- feature_count[nonzero_idx] + 1L
      best_model <- list(test_auc = as.numeric(roc_test$auc),
                         model = cv_fit, lambda = lambda,
                         train_idx = idx, test_idx = setdiff(seq_len(n_total), idx))
    }

    if (iterative > 1 && ok) {
      success_n <- success_n + 1L
      feature_count[nonzero_idx] <- feature_count[nonzero_idx] + 1L
      if (as.numeric(roc_test$auc) > best_model$test_auc)
        best_model <- list(test_auc = as.numeric(roc_test$auc),
                           model = cv_fit, lambda = lambda,
                           train_idx = idx, test_idx = setdiff(seq_len(n_total), idx))
    }
  }
  cli::cli_progress_done()

  # frequency table and final performance
  feature_freq <- dplyr::tibble(Feature = names(feature_count),
                                count = as.integer(feature_count)) %>% dplyr::arrange(dplyr::desc(count))
  if (!is.null(best_model$model)) {
    train_df_best <- base_df[best_model$train_idx, ]; test_df_best <- base_df[best_model$test_idx, ]
    x_train_best <- train_df_best %>% dplyr::select(-PHENOTYPE) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>% as.matrix()
    y_train_best <- train_df_best$PHENOTYPE
    x_test_best  <- test_df_best  %>% dplyr::select(-PHENOTYPE) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>% as.matrix()
    y_test_best  <- test_df_best$PHENOTYPE

    s_train <- drop(stats::predict(best_model$model, newx = x_train_best, s = best_model$lambda, type = score_type))
    s_test  <- drop(stats::predict(best_model$model, newx = x_test_best,  s = best_model$lambda, type = score_type))
    perf_train <- prs.roc_hist_density(s_train, y_train_best, ctrl_case_level = c(0, 1))
    perf_test  <- prs.roc_hist_density(s_test,  y_test_best,  ctrl_case_level = c(0, 1))
  } else { perf_train <- NULL; perf_test <- NULL }

  leo.basic::leo_log("iLasso done: success={success_n}/{attempts}, best_test_auc={round(max(auc_history[auc_history$pass,]$test_auc, na.rm = T), 3)}", level = "success")

  list(model = best_model$model, lambda = best_model$lambda,
       perf_train = perf_train, perf_test = perf_test,
       success_n = success_n, attempts = attempts,
       auc_history = auc_history, feature_rank = feature_freq)
}

#' @rdname lasso_prs
#' @export
lasso_prs_target <- function(a1_matrix, model, lambda, score_type = "link"){
  df <- a1_matrix %>% dplyr::mutate(PHENOTYPE = ifelse(PHENOTYPE == 1, 0L, 1L)) %>% dplyr::select(-c(FID, IID, PAT, MAT, SEX))
  x_test <- df %>% dplyr::select(-PHENOTYPE) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>% as.matrix()
  y_test <- df$PHENOTYPE

  pred_score <- drop(stats::predict(model, newx = x_test, s = lambda, type = score_type))
  perf <- prs.roc_hist_density(pred_score, y_test, ctrl_case_level = c(0, 1))

  pred_df <- data.frame(FID = a1_matrix$FID, PHENOTYPE = y_test, PRED_SCORE = pred_score)
  return(list(pred_df = pred_df, perf = perf,
              target_x = x_test, target_y = y_test))
}
#' @rdname lasso_prs
#' @export
lasso_prs_rank <- function(model, rank, auc_history){
  return_list <- list()
  # rank plot
  rank_plot <- rank %>%
    arrange(count) %>%
    mutate(Feature = factor(Feature, levels = Feature),
           Occurence = count) %>%
    ggplot(aes(x = Occurence, y = Feature)) +
    geom_bar(stat = "identity", fill = "maroon") +
    labs(title = "iLasso Feature Occurence", x = NULL, y = NULL) +
    theme_bw() +
    # scale_x_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
    theme(plot.title.position = "panel", # or "plot"
          plot.title = element_text(hjust = 0, size = 20, face = "bold"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 13, color = "black"),
          axis.line = element_line(linewidth = .5, color = "black"),
          legend.position = c(0.85, 0.7), # Place legend inside the plot at (0.8, 0.8)
          legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
          legend.key = element_blank(), # Remove border around legend keys
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank() )
  return_list$rank_plot <- rank_plot

  # train value distribution
  plot_p_u_test <- ggplot(auc_history, aes(x = iteration, y = -log10(p_u_test)))+
    stat_binhex(alpha=1, bins = 20)+
    scale_fill_gradient(low = "lightgray",high = "maroon")+
    xlab("Iteration")+ylab('-Log10 P Value (Test Set)')+
    geom_rug(position = "jitter", size = 0.1, color = "black", alpha = 0.75)+
    theme_classic() +
    theme(plot.title.position = "panel", # or "plot"
          plot.title = element_text(hjust = 0, size = 20, face = "bold"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 13, color = "black"),
          axis.line = element_line(linewidth = .5, color = "black"),
          legend.position = "right", # Place legend inside the plot at (0.8, 0.8)
          legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
          legend.key = element_blank(), # Remove border around legend keys
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank() )
  plot_p_delong <- ggplot(auc_history, aes(x = iteration, y = p_delong))+
    stat_binhex(alpha = 1, bins = 20)+
    scale_fill_gradient(low = "lightgray",high = "maroon")+
    xlab("Iteration")+ylab('-Log10 P Value (Test Set)')+
    geom_rug(position = "jitter", size = 0.1, color = "black", alpha = 0.75)+
    theme_classic() +
    theme(plot.title.position = "panel", # or "plot"
          plot.title = element_text(hjust = 0, size = 20, face = "bold"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 13, color = "black"),
          axis.line = element_line(linewidth = .5, color = "black"),
          legend.position = "right", # Place legend inside the plot at (0.8, 0.8)
          legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
          legend.key = element_blank(), # Remove border around legend keys
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank() )
  plot_n_feat <- ggplot(auc_history, aes(x = iteration, y = n_feat))+
    stat_binhex(alpha = 1, bins = 20)+
    scale_fill_gradient(low = "lightgray",high = "maroon")+
    xlab("Iteration")+ylab('Feature Number')+
    geom_rug(position = "jitter", size = 0.1, color = "black", alpha = 0.75)+
    theme_classic() +
    theme(plot.title.position = "panel", # or "plot"
          plot.title = element_text(hjust = 0, size = 20, face = "bold"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 13, color = "black"),
          axis.line = element_line(linewidth = .5, color = "black"),
          legend.position = "right", # Place legend inside the plot at (0.8, 0.8)
          legend.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
          legend.key = element_blank(), # Remove border around legend keys
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank() )
  return_list$plot_p_u_test <- plot_p_u_test
  return_list$plot_p_delong <- plot_p_delong
  return_list$plot_n_feat <- plot_n_feat
  return(return_list)
}

# plink ----
plink_get_a1_matrix <- function(bfile, summary, output_dir, output_name,
                                snp_col="SNP.meta", a1_col="A1",
                                plink_bin = plinkbinr::get_plink_exe()) {
  requireNamespace("data.table"); requireNamespace("dplyr")
  if (!file.exists(paste0(bfile, ".bed"))) stop("Base .bed file not found.")
  if (!file.exists(paste0(bfile, ".bim"))) stop("Base .bim file not found.")
  if (!file.exists(paste0(bfile, ".fam"))) stop("Base .fam file not found.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.character(summary)) {
    summary_file <- fread(summary) %>% dplyr::transmute(SNP=.data[[snp_col]], A1=.data[[a1_col]])
  } else {
    summary_file <- summary %>% dplyr::transmute(SNP=.data[[snp_col]], A1=.data[[a1_col]])
  }
  # output_name <- paste0(output_name, "_", format(Sys.time(), "%Y%m%d"))
  snp_to_extract <- file.path(output_dir, paste0(output_name, ".snp_to_extract.txt"))
  snp_a1_txt  <- file.path(output_dir, paste0(output_name, ".snp_a1.txt"))
  summary_file %>% dplyr::select(SNP) %>% fwrite(snp_to_extract, sep="\t", col.names=F)
  summary_file %>% dplyr::select(SNP, A1) %>% fwrite(snp_a1_txt, sep="\t", col.names=F)
  leo.basic::leo_log("-----\nSNP list for extraction saved to {.path {snp_to_extract}}\nSNP-A1 saved to {.path {snp_a1_txt}}")

  cmd <- paste(
    plink_bin,
    "--bfile", bfile,
    "--extract", snp_to_extract,
    "--recode A include-alt",
    "--recode-allele", snp_a1_txt,
    "--out", file.path(output_dir, output_name)
  )
  system(cmd)
  leo.basic::leo_log("Genotype matrix extracted: {.path {file.path(output_dir, output_name)}.raw}", level="success")
  fread(paste0(file.path(output_dir, output_name),".raw"))
}

plink_prs <- function(bfile, snp_file, output_prefix, plink_bin = plinkbinr::get_plink_exe()) {
  # This function calculate the PRS using PLINK
  if (!file.exists(paste0(bfile, ".bed"))) stop("Base .bed file not found.")
  if (!file.exists(paste0(bfile, ".bim"))) stop("Base .bim file not found.")
  if (!file.exists(paste0(bfile, ".fam"))) stop("Base .fam file not found.")
  if (!file.exists(snp_file)) stop("SNP file not found.")
  cmd <- paste(
    plink_bin,
    "--bfile", bfile,
    "--score", snp_file, "1", "2", "3", "header", "sum",
    "--out", output_prefix
  )
  system(cmd)
  leo.basic::leo_log("PRS calculation completed. Output saved to {.path {output_prefix}.profile}", level = "success")
}

plink_prs_target <- function(bfile, snp_file, output_prefix, plink_bin = plinkbinr::get_plink_exe()){
  # This function evaluate
  stage3_plink <- plink_prs(bfile = "./data/zuo/bed/SNP_for_PRS/stage3/VKH-ASA-141s-1012s-ALL-PRS",
                            snp_file = snp_file, output_prefix = output_prefix, plink_bin = plink_bin)
  stage3_plink_p <- paste(file.path(output_dir, "stage3_plink"), "profile", sep = ".")
  stage3_prs_plink <- fread(stage3_plink_p)
  stage3_prs_plink <- stage3_prs_plink %>% mutate(PHENO = ifelse(PHENO == 1, 0L, 1L))
}

# plink clump ----
#' HLA-aware PLINK clumping (Yet implemented)
#'
#' Simple wrapper:
#'   "none"       -> drop HLA (chr6:25–34Mb), clump non-HLA only
#'   "all"        -> keep all HLA (no clump), clump non-HLA, then union
#'   "just_clump" -> clump all variants together
#'
#' Input: `{output}.clump.txt` with columns "SNP","P".
#' Defaults follow common PRS C+T: --clump-p1 1, --clump-p2 1, --clump-r2 0.1, --clump-kb 250.
#'
#' @param bfile PLINK bfile for LD reference.
#' @param hla_keep One of c("all","none","just_clump").
#' @param output Output prefix; expects `{output}.clump.txt`.
#' @param plink_bin Path to PLINK.
#'
#' @return Character vector of kept SNPs; also writes `{output}.kept.snps.txt`.
#' @export
#' @importFrom data.table fread fwrite
#' @importFrom dplyr left_join select transmute
#' @importFrom stats setNames
plink_clump_hla_aware <- function(bfile, hla_keep = c("all","none","just_clump"),
                                  output, plink_bin = plinkbinr::get_plink_exe()) {
  # basic checks
  if (!file.exists(paste0(bfile, ".bed"))) stop("Base .bed file not found.")
  if (!file.exists(paste0(bfile, ".bim"))) stop("Base .bim file not found.")
  if (!file.exists(paste0(bfile, ".fam"))) stop("Base .fam file not found.")
  hla_keep <- match.arg(hla_keep)

  in_file <- paste0(output, ".clump.txt"); out_pref <- paste0(output, ".run")
  if (!file.exists(in_file)) stop("Input file not found: {output}.clump.txt (need columns: SNP, P)")

  # annotate CHR/BP to mark HLA region
  bim <- data.table::fread(paste0(bfile, ".bim"), header = F, col.names = c("CHR","SNP","CM","BP","A1","A2"))
  in_df <- data.table::fread(in_file)
  if (!all(c("SNP","P") %in% colnames(in_df))) stop("Input must have columns: SNP, P")
  in_df <- dplyr::left_join(in_df, bim %>% dplyr::select(CHR, SNP, BP), by = "SNP")
  in_hla <- with(in_df, CHR == 6 & BP >= 25e6 & BP <= 34e6)
  df_hla <- in_df[in_hla, , drop = F]; df_non <- in_df[!in_hla, , drop = F]

  # helper: run PLINK clump on 2-col file (SNP,P)
  run_clump <- function(df, tag) {
    if (nrow(df) == 0) return(character(0))
    tag_in <- paste0(out_pref, ".", tag, ".in.txt"); tag_out <- paste0(out_pref, ".", tag)
    data.table::fwrite(df %>% dplyr::transmute(SNP, P = as.numeric(P)), tag_in, sep = "\t", quote = F, col.names = T)
    cmd <- paste(plink_bin, "--bfile", bfile, "--clump", tag_in,
                 "--clump-snp-field", "SNP", "--clump-p-field", "P",
                 "--clump-p1", 1, "--clump-p2", 1, "--clump-r2", 0.1, "--clump-kb", 250, "--out", tag_out)
    base::system(cmd)
    clumped <- paste0(tag_out, ".clumped")
    if (!file.exists(clumped)) { leo.basic::leo_log("No .clumped for {tag}; return empty.", level = "warning"); return(character(0)) }
    tbl <- tryCatch(data.table::fread(clumped, fill = T), error = function(e) NULL)
    if (is.null(tbl) || !"SNP" %in% colnames(tbl)) { leo.basic::leo_log("Parse .clumped failed for {tag}.", level = "danger"); return(character(0)) }
    unique(as.character(tbl$SNP))
  }

  # modes
  if (hla_keep == "none") {
    keep <- run_clump(df_non %>% dplyr::select(SNP, P), "nonHLA")
  } else if (hla_keep == "all") {
    keep_non <- run_clump(df_non %>% dplyr::select(SNP, P), "nonHLA")
    keep <- unique(c(keep_non, df_hla$SNP))
  } else {
    keep <- run_clump(in_df %>% dplyr::select(SNP, P), "all")
  }

  data.table::fwrite(data.frame(SNP = keep), paste0(output, ".kept.snps.txt"), sep = "\t", quote = F, col.names = T)
  leo.basic::leo_log("Clump done ({hla_keep}): input={nrow(in_df)}, kept={length(keep)}; saved -> {output}.kept.snps.txt", level = "success")
  return(keep)
}
# vis ----
prs.roc_hist_density <- function(prs_score, phenotype,
                                 ctrl_case_level = c(1,2),
                                 bins = 100,
                                 palette2 = c("#4C78A8", "#F58518")) {
  stopifnot(length(prs_score) == length(phenotype))

  # 1. Clean phenotype code --> 0/1
  control_code <- ctrl_case_level[1]
  case_code    <- ctrl_case_level[2]
  phenotype_binary <- ifelse(phenotype == case_code, 1L,
                             ifelse(phenotype == control_code, 0L, NA_integer_))

  valid_index <- complete.cases(prs_score, phenotype_binary)
  prs_score   <- prs_score[valid_index]
  phenotype_binary <- phenotype_binary[valid_index]

  # 2. ROC
  roc_obj <- pROC::roc(phenotype_binary, prs_score)
  cli::cli_alert_info("AUC = {roc_obj$auc}")
  roc_df <- data.frame(TPR = rev(roc_obj$sensitivities),
                       FPR = rev(1 - roc_obj$specificities))
  roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = "black", linewidth = .7) +
    geom_ribbon(aes(ymin = 0, ymax = TPR), fill = ggsci::pal_npg("nrc")(1), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linewidth = .75, linetype = "dashed", color = "black") +
    annotate("text", x = 0.75, y = 0.1, size = 6, color = "black", # fontface = "bold",
             label = paste("AUC =", sprintf("%.3f", roc_obj$auc))) +
    labs(title = "ROC Curve",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_classic() +
    # scale_x_continuous(expand = c(0, 0.01), labels = function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", format(x, digits = 2)))) +
    # scale_y_continuous(expand = c(0, 0.01), labels = function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", format(x, digits = 2)))) +
    theme(
      plot.title.position = "panel", # or "plot"
      plot.title = element_text(hjust = 0, size = 20, face = "bold"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.line = element_line(linewidth = .5, color = "black"),
      legend.position = c(0.85, 0.7), # Place legend inside the plot at (0.8, 0.8)
      legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
      legend.key = element_blank(), # Remove border around legend keys
      legend.text = element_text(size = 12)
    )

  # hist plot
  res <- data.frame(Phenotype = phenotype_binary, PRS = prs_score)
  medians <- res %>%
    group_by(Phenotype) %>%
    summarise(median = median(PRS, na.rm = TRUE), .groups = "drop")
  hist <- ggplot(res, aes(x = prs_score, fill = factor(Phenotype))) +
    geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
    geom_vline(data = medians, aes(xintercept = median, color = factor(Phenotype)),
               linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = c(palette2), labels = c("Control", "Case")) +
    scale_color_manual(values = c(palette2), guide = "none") +
    labs(title = "PRS Distribution", x = "PRS", y = "Count", fill = "Phenotype") +
    theme_classic() +
    # scale_x_continuous(expand = c(0, 0.01)) +
    expand_limits(y = 0) +
    theme(
      # panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0, size = 20, face = "bold"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.line = element_line(linewidth = .5, color = "black"),
      legend.position = c(0.85, 0.7), # Place legend inside the plot at (0.8, 0.8)
      legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
      legend.background = element_rect(fill = "gray90"), # Add lightgray background to legend
      legend.key = element_blank(), # Remove border around legend keys
      legend.text = element_text(size = 12)
    )

  # density plot
  dens_plot <- ggplot(res, aes(x = PRS, fill = factor(Phenotype))) +
    geom_density(alpha = .5, size = .7) +
    scale_fill_manual(values = palette2, labels = c("Control", "Case")) +
    labs(title = "PRS Distribution", x = "PRS Score", y = "Density",  fill = "Phenotype") +
    theme_classic() +
    # scale_x_continuous(expand = c(0, 0.001)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(
      # panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0, size = 20, face = "bold"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.line = element_line(linewidth = .5, color = "black"),
      legend.position = c(0.85, 0.7), # Place legend inside the plot at (0.8, 0.8)
      legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
      legend.background = element_rect(fill = "gray90"), # Add lightgray background to legend
      legend.key = element_blank(), # Remove border around legend keys
      legend.text = element_text(size = 12)
    )
  return(list(roc_value = roc_obj$auc,
              roc_plot  = roc_plot,
              hist_plot = hist,
              dens_plot = dens_plot))
}

gglasso <- function(no_cv_lasso, cv_lasso, pos = "right", guide_cols = 1) {
  # no_cv_lasso: glmnet fit returned by glmnet::glmnet()
  # cv_lasso   : cv.glmnet object (provides lambda.min / lambda.1se)
  # pos        : legend position ("none","left","right","bottom","top")
  # guide_cols : number of legend columns

  cm <- as.matrix(stats::coef(no_cv_lasso))                  # coef matrix: rows=features, cols=path steps
  rn <- rownames(cm)
  if ("(Intercept)" %in% rn) cm <- cm[setdiff(rn, "(Intercept)"), , drop = FALSE]

  colnames(cm) <- paste0("s", seq_len(ncol(cm)))
  rn <- rownames(cm)

  df <- tibble::as_tibble(cm, .name_repair = "check_unique") %>% dplyr::mutate(coef = rn)
  long <- tidyr::pivot_longer(df, -coef, names_to = "s", values_to = "value") %>%
    dplyr::mutate(idx = as.integer(sub("^s", "", s)),
                  lambda = no_cv_lasso$lambda[idx],
                  coef = factor(coef))

  # optional L1 norm if needed:
  # l1_norm <- colSums(abs(cm[setdiff(rownames(cm), "(Intercept)"), , drop = FALSE]))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = log(lambda), y = value, color = coef)) +
    ggplot2::geom_vline(xintercept = log(cv_lasso$lambda.min), size = 0.8, color = "grey60", alpha = 0.8, linetype = 2) +
    ggplot2::geom_vline(xintercept = log(cv_lasso$lambda.1se), size = 0.8, color = "grey60", alpha = 0.8, linetype = 2) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::xlab("Lambda (log scale)") + ggplot2::ylab("Coefficients") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 12, color = "black"),
                   axis.title = ggplot2::element_text(size = 15, color = "black"),
                   axis.text.x = ggplot2::element_text(margin = grid::unit(c(0.25, 0.25, 0.25, 0.25), "cm")),
                   axis.text.y = ggplot2::element_text(margin = grid::unit(c(0.25, 0.25, 0.25, 0.875), "cm")),
                   axis.title.y = ggplot2::element_text(vjust = -1),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 12, color = "black"),
                   legend.position = pos,
                   plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = guide_cols))
  return(p)
}

ggcvlasso <- function(cv_lasso, pos = "right", guide_cols = 1) {
  # cv_lasso   : cv.glmnet object
  # pos        : legend position
  # guide_cols : number of legend columns

  xx <- data.frame(lambda = cv_lasso$lambda, cvm = cv_lasso$cvm, cvsd = cv_lasso$cvsd,
                   cvup = cv_lasso$cvup, cvlo = cv_lasso$cvlo, nzero = cv_lasso$nzero)
  xx$ll <- log(xx$lambda)
  xx$NZERO <- ifelse(xx$nzero == 1, "1 biomarker", paste0(xx$nzero, " biomarkers"))
  xx$NZERO <- factor(xx$NZERO, levels = rev(unique(xx$NZERO)))
  ylab_txt <- if (!is.null(cv_lasso$name)) cv_lasso$name else "CV metric"

  p <- ggplot2::ggplot(xx, ggplot2::aes(ll, cvm, color = NZERO)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = cvlo, ymax = cvup), width = 0.05, linewidth = 1) +
    ggplot2::geom_vline(xintercept = log(cv_lasso$lambda.min), size = 0.8, color = "grey60", alpha = 0.8, linetype = 2) +
    ggplot2::geom_vline(xintercept = log(cv_lasso$lambda.1se), size = 0.8, color = "grey60", alpha = 0.8, linetype = 2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::xlab("Lambda (log scale)") + ggplot2::ylab(ylab_txt) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0.02)) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 12, color = "black"),
                   axis.title = ggplot2::element_text(size = 15, color = "black"),
                   axis.text.x = ggplot2::element_text(margin = grid::unit(c(0.25, 0.25, 0.25, 0.25), "cm")),
                   axis.text.y = ggplot2::element_text(margin = grid::unit(c(0.25, 0.25, 0.25, 1), "cm")),
                   axis.title.y = ggplot2::element_text(vjust = -1),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 12, color = "black"),
                   legend.position = pos,
                   plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = guide_cols))
  return(p)
}

# helper ----
.catboost_install_message <- function(context = NULL){
  context_txt <- if (is.null(context)) "" else paste0(context, "\n")
  paste0(context_txt,
         "Package 'catboost' is required.\n",
         "Install (binary, recommended as of 2026-02):\n",
         "  pak::pak('url::https://github.com/catboost/catboost/releases/download/v1.2.8/catboost-R-darwin-universal2-1.2.8.tgz')\n",
         "See: https://catboost.ai/docs/en/installation/r-installation-binary-installation")
}

.check_catboost <- function(context = NULL){
  if (!requireNamespace("catboost", quietly = TRUE)) stop(.catboost_install_message(context = context), call. = FALSE)
}

.check_caret <- function(context = NULL){
  if (!requireNamespace("caret", quietly = TRUE)) stop(paste0("Package 'caret' is required", if (is.null(context)) "." else paste0(" (", context, ")."), " Install with install.packages('caret')."), call. = FALSE)
}
