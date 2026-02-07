# find_proxy

`find_proxy` finds the proxy snp for the miss iv

## Usage

``` r
find_proxy(
  miss_iv,
  miss_snp,
  outcome_snp,
  proxy_file = NULL,
  proxy_output_path = NULL,
  pop = "EUR",
  gb = "grch38",
  token = ""
)
```

## Arguments

- miss_iv:

  iv datafram in tsmr exposure format, which can not locate snp in the
  outcome

- miss_snp:

  snp in the miss_iv; can be inferred using miss_iv\$SNP

- outcome_snp:

  a str vector containing all snp in the outcome; this NOT the entire
  outcome gwas summary statistics

- proxy_file:

  pre-calculated proxy file path (full); do provide this if the proxy
  file is already generated !!!

- proxy_output_path:

  a full path to save the proxy file when using ldlink

- pop:

  reference panel from 1kg (LDlinkR param)

- gb:

  genome build (LDlinkR param)

- token:

  token of `LDlinkR`

## Value

a updated missiv with `proxy.snp` `proxy.effect.allele`
`proxy.other.allele` `r2` col

## Examples

``` r
# This function can be used when many iv can not locate corresponding snp in the outcome in tsmr analysis
miss_iv <- iv[!iv$SNP %in% dat_h$SNP,] # iv is estracted iv via tsmr package;dat_h is a standard output of harmonise_data()
#> Error: object 'iv' not found
miss_snp <- miss_iv$SNP
#> Error: object 'miss_iv' not found
outcome_snp <- iri_nc$SNP
#> Error: object 'iri_nc' not found
proxy_output_path <- "Full path to where you wanna store the LDlinkR output"
 proxy_iv <- find_proxy(miss_iv, miss_snp, outcome_snp,
             proxy_file = "/Users/leoarrow/project/iridocyclitis/output/tsmr//combined_query_snp_list_grch38.txt",
             proxy_output_path = NULL)
#> Error in find_proxy(miss_iv, miss_snp, outcome_snp, proxy_file = "/Users/leoarrow/project/iridocyclitis/output/tsmr//combined_query_snp_list_grch38.txt",     proxy_output_path = NULL): could not find function "find_proxy"
 # bak
 proxy_iv$target.snp <- proxy_iv$SNP # target snp
#> Error: object 'proxy_iv' not found
 proxy_iv$target.A1 <- proxy_iv$effect_allele.exposure
#> Error: object 'proxy_iv' not found
 proxy_iv$target.A2 <- proxy_iv$other_allele.exposure
#> Error: object 'proxy_iv' not found
 # replace for tsmr
 proxy_iv$SNP <- proxy_iv$proxy.snp
#> Error: object 'proxy_iv' not found
 proxy_iv$effect_allele.exposure <- proxy_iv$proxy.A1
#> Error: object 'proxy_iv' not found
 proxy_iv$other_allele.exposure <- proxy_iv$proxy.A2
#> Error: object 'proxy_iv' not found
 iv_f <- bind_rows(non_miss_iv, proxy_iv) # f for final
#> Error: object 'non_miss_iv' not found
 dat_h_proxy <- harmonise_data(iv_f, out_nc_proxy)
#> Error in harmonise_data(iv_f, out_nc_proxy): could not find function "harmonise_data"
 mr(dat_h_proxy) # nailed it!
#> Error in mr(dat_h_proxy): could not find function "mr"
```
