# Iterative Lasso PRS (iLasso) utilities

Train, apply and visualize iterative Lasso-based PRS models with success
gating.

## Usage

``` r
lasso_prs(
  a1_matrix,
  divide_ratio = 0.7,
  iterative = 100,
  nfolds = 10,
  score_type = "link",
  seed = 725
)

lasso_prs_target(a1_matrix, model, lambda, score_type = "link")

lasso_prs_rank(model, rank, auc_history)
```

## Arguments

- a1_matrix:

  data.frame; PLINK A1 matrix.

- divide_ratio:

  numeric; train fraction per iteration (default 0.7).

- iterative:

  integer; number of iterations (default 100) - If set to 1, then
  regular lasso is performed.

- nfolds:

  integer; CV folds for glmnet (default 10).

- score_type:

  character; "link" (linear score, default) or "response" (probability).

- seed:

  integer; base random seed (default 725).

- model:

  Trained glmnet model (for target/rank functions).

- lambda:

  Best lambda from CV.

- rank:

  Feature frequency table.

- auc_history:

  Tibble of iteration results.

## Value

Each function returns a list:

- `lasso_prs`: model, lambda, perf_train, perf_test, success_n,
  attempts, auc_history, feature_rank

- `lasso_prs_target`: pred_df, perf, target_x, target_y

- `lasso_prs_rank`: plots for feature frequency and diagnostics

## See also

[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
[`roc`](https://rdrr.io/pkg/pROC/man/roc.html)
