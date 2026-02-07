# HLA-aware PLINK clumping (Yet implemented)

Simple wrapper: "none" -\> drop HLA (chr6:25–34Mb), clump non-HLA only
"all" -\> keep all HLA (no clump), clump non-HLA, then union
"just_clump" -\> clump all variants together

## Usage

``` r
plink_clump_hla_aware(
  bfile,
  hla_keep = c("all", "none", "just_clump"),
  output,
  plink_bin = plinkbinr::get_plink_exe()
)
```

## Arguments

- bfile:

  PLINK bfile for LD reference.

- hla_keep:

  One of c("all","none","just_clump").

- output:

  Output prefix; expects `{output}.clump.txt`.

- plink_bin:

  Path to PLINK.

## Value

Character vector of kept SNPs; also writes `{output}.kept.snps.txt`.

## Details

Input: `{output}.clump.txt` with columns "SNP","P". Defaults follow
common PRS C+T: –clump-p1 1, –clump-p2 1, –clump-r2 0.1, –clump-kb 250.
