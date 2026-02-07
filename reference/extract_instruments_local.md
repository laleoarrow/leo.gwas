# Extract instruments locally for MR Analysis

Filters SNPs by p-value threshold and performs LD clumping with flexible
column mapping.

## Usage

``` r
extract_instruments_local(
  dat,
  p = 5e-08,
  pop = NULL,
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
  plink_bin = NULL
)
```

## Arguments

- dat:

  Dataframe containing GWAS summary statistics; change it to dataframe
  if it is data.table.

- p:

  P-value cutoff (default = 5e-8)

- pop:

  Super-population code for LD reference (default = NULL, as we
  recommend using loca bfile)

- phenotype_col:

  Column name for phenotype (default = "Phenotype")

- snp_col:

  Column name for SNP IDs (default = "SNP")

- chr_col:

  Column name for chromosome (default = "CHR")

- pos_col:

  Column name for position (default = "POS")

- effect_allele_col:

  Column name for effect allele (default = "A1")

- other_allele_col:

  Column name for non-effect allele (default = "A2")

- beta_col:

  Column name for effect size (default = "BETA")

- se_col:

  Column name for standard error (default = "SE")

- pval_col:

  Column name for p-values (default = "P")

- eaf_col:

  Column name for effect allele frequency (default = "EAF")

- N:

  Column name for sample size (default = "Neff", set NULL to exclude)

- bfile:

  Path to PLINK binary reference panel (default = NULL)

- plink_bin:

  Path to PLINK executable (default = NULL)

## Value

Formatted exposure data ready for MR analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# Custom column names example
extract_instruments_local(
  dat = gwas_data,
  phenotype_col = "Trait",
  snp_col = "rsID",
  chr_col = "Chromosome",
  pos_col = "Position"
)
} # }
```
