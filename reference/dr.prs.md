# DuoRank PRS (Dr.PRS)

Dr.PRS is designed to rank the importance for PRS inputs and futher
optimazation. It takes two none-overlap stage plink file (bfile) for
machine learning modeling:

- lasso (to capture linear relationship)

- catboost (to capture non-linear relationship) It also calculated PRS
  with plink using traditional additive model.

## Usage

``` r
dr.prs(
  stage1_bfile = "./data/zuo/bed/SNP_for_PRS/VKH-zhonghua-for_prs_snplist",
  stage2_bfile = "./data/zuo/bed/SNP_for_PRS/VKH-ASA-for_prs_snplist",
  summary_file = "./output/part1-gwas/prs/snp_for_prs.txt",
  output_dir = "./output/part1-gwas/prs/b1t2",
  method = c(1, 2, 3),
  dr_optimazation = TRUE,
  snp_col = "SNP.meta",
  a1_col = "A1",
  weight_col = "OR_Final",
  weight_type = "OR",
  divide_ratio = 0.7,
  iLasso_iteration = 1000,
  nfolds = 10,
  plink_bin = plinkbinr::get_plink_exe(),
  clump = T,
  clump_include_hla = F,
  clump_param = NULL,
  seed = 725
)

combine_rank(rank1, rank2, auc1 = NULL, auc2 = NULL)
```

## Arguments

- stage1_bfile:

  Path to the stage1 PLINK binary files. Normally it is the one
  generates summary data. If only 1 source of individual data is
  available, you can split it into 2 non-overlap datasets.

- stage2_bfile:

  Path to the stage2 PLINK binary files, i.e., the target dataset for
  PRS calculation.

- summary_file:

  Path to the summary statistics file or a data frame containing SNPs,
  alleles, and weights. Note this summary only contains refined SNP for
  PRS calculation.

- output_dir:

  Directory to save output files.

- method:

  Integer vector indicating which methods to use: 1=PLINK, 2=CatBoost,
  3=iLasso. Default is all three methods.

- dr_optimazation:

  Logical, whether to perform greedy search to optimize PRS based on
  combined rank from CatBoost and iLasso (default: TRUE).

- snp_col:

  Column name in the summary statistics for SNP IDs (default:
  "SNP.meta").

- a1_col:

  Column name in the summary statistics for effect allele (default:
  "A1").

- weight_col:

  Column name in the summary statistics for weights (default:
  "OR_Final").

- weight_type:

  Type of weights, either "OR" (default) or "Beta". If "OR", it will
  calculate Beta values.

- divide_ratio:

  Proportion of stage1 data in iLasso (default: 0.7).

- iLasso_iteration:

  Number of iterations for iLasso model (default: 1000).

- nfolds:

  Number of folds for cross-validation (default: 10).

- plink_bin:

  Path to the PLINK binary executable (default:
  `plinkbinr::get_plink_exe()`).

- clump:

  Place holder for future clumping function.

- clump_include_hla:

  Place holder for future clumping function.

- clump_param:

  Place holder for future clumping function.

- seed:

  Random seed (Default: 725).

- rank1:

  Data frame containing rank information for the first stage.

- rank2:

  Data frame containing rank information for the second stage.

- auc1:

  Numeric. AUC for the first stage.

- auc2:

  Numeric. AUC for the second stage.

## Value

A list with:

- lasso_importance, catboost_importance, combined_rank

- greedy_auc_curve (data.frame: kept_snps, auc)

- final_subset (character vector of SNPs)

- plink_prs_file (path), logs
