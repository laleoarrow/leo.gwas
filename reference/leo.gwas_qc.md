# GWAS summary QC pipeline (chip + imputed)

Perform QC on GWAS summary stats by merging imputed and genotyped data,
removing duplicates, checking consistency, matching to a reference
panel, filtering by DAF/F_U cutoffs, and saving results.

## Usage

``` r
leo.gwas_qc(
  summary_x2_p,
  summary_x2_chip_p,
  ref_p = "/Users/leoarrow/Project/ref/1kg_maf/zuo_ref/1KG-EAS-EAF-chrposa2a1.txt.gz",
  save_dir = "~/Project/BD2025/data/qc/sex_split",
  save_name_prefix = "bd-ASA-41",
  DAF_cutoff = 0.2,
  F_U_cutoff = 0.01
)
```

## Arguments

- summary_x2_p:

  Path to imputed summary file (.assoc/.gz).

- summary_x2_chip_p:

  Path to genotyped (chip) summary file.

- ref_p:

  Path to reference panel or loaded df

- save_dir:

  Output directory (default: `"~/Project/BD2025/data/qc/sex_split"`).

- save_name_prefix:

  Prefix for output file names.

- DAF_cutoff:

  Numeric. Max \|DAF\| allowed (default 0.2).

- F_U_cutoff:

  Numeric. Min F_U required (default 0.01).

## Value

return(summary_qc); writes QC results to `save_dir`.

## Details

Major steps:

- Read and standardize imputed & chip data

- Remove duplicates (keep smallest P for genotyped, drop multiple
  imputed)

- Merge and label GI (imputed/genotyped)

- Match with reference panel and compute DAF

- Filter by cutoffs and save full / P\<1e-6 subsets and output

Logs at each step with `leo_log()` for dimension, duplication, NA, etc.

## Examples

``` r
if (FALSE) { # \dontrun{
leo.gwas_qc("imp.assoc.gz", "chip.assoc")
} # }
```
