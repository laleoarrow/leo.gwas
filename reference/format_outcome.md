# format outcome data

format outcome data

## Usage

``` r
format_outcome(dat, snp = iv$SNP, N = "Neff")
```

## Arguments

- dat:

  a dataframe for outcome with SNP, CHR, POS, A1, A2, EAF, BETA, SE, P,
  Phenotype, samplesize columns

- snp:

  a str vector out of iv\$SNP

- N:

  N column name for sample size (effective or observed)

## Value

a tsmr format outcome dataframe
