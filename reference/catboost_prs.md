# CatBoost PRS utilities

Train, apply, and rank CatBoost polygenic risk score (PRS) models.

## Usage

``` r
catboost_prs(a1_matrix, divide = F, divide_ratio = 0.5)

catboost_prs_target(a1_matrix, model)

catboost_prs_rank(
  model,
  pool,
  pool_df,
  types = c("FeatureImportance", "ShapValues", "Interaction"),
  top_k = NULL
)
```

## Arguments

- a1_matrix:

  A1 matrix (with PHENOTYPE and ID columns).

- divide:

  Logical; split into train/val sets. Default FALSE.

- divide_ratio:

  Numeric; train proportion if divide=TRUE. Default 0.5.

- model:

  Trained CatBoost model.

- pool:

  CatBoost pool.

- pool_df:

  Data used to build `pool`.

- types:

  Character; any of "FeatureImportance","ShapValues","Interaction".

- top_k:

  Integer; top features to keep. Default NULL.

## Value

Each function returns a list; contents depend on task (training, target
prediction, ranking).

## Details

These functions require installed catboost. Training with
`divide = TRUE` also requires caret for stratified data splitting.

## See also

`catboost::catboost.train()`,
`catboost::catboost.get_feature_importance()`
