# mrlap_one_pair

`mrlap_one_pair` perform the MR-lap for one pair of exposure and outcome
[MR-lap](https://github.com/n-mounier/MRlap?tab=readme-ov-file)

## Usage

``` r
mrlap_one_pair(
  exposure_data,
  outcome_data,
  exposure_name,
  outcome_name,
  ld_path,
  hm3_path,
  log = F,
  log_path = "."
)
```

## Arguments

- exposure_data:

  exposure_data

- outcome_data:

  outcome_data

- exposure_name:

  exposure_name

- outcome_name:

  outcome_name

- ld_path:

  ldsc required file

- hm3_path:

  ldsc required file; tutorial use nomhc version list

- log:

  logi, if T, would save ldsc log to log_path, defaut F

- log_path:

  path to where you wanna store the ldsc log
