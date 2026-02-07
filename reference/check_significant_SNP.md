# locate the significant SNP for conditional analysis

locate the significant SNP for conditional analysis

## Usage

``` r
check_significant_SNP(x, environment.list, significance_level = 5e-08)
```

## Arguments

- x:

  data.frame of the SNP information

- environment.list:

  tibble of the environment list (see the example for usage)

- significance_level:

  numeric, significance level, default is 5e-8

## Value

message that informs the user of the independant SNP for the following
conditional analysis

## HLA data analysis

## Examples

``` r
con_dir <- "/Users/leoarrow/project/VKH2024/data/zuo/con_su" # specify the directory to store the HLA original data and subsequent conditional analysis results data
files <- list.files(con_dir,full.names = T) %>% as.vector(); files # update it each time
#> Error in list.files(con_dir, full.names = T) %>% as.vector(): could not find function "%>%"
x1 <- fread(files[1]) %>% arrange(desc(CHISQ)); head(x1) # ! for the first one, just read the data and sort it by CHISQ/P value
#> Error in fread(files[1]) %>% arrange(desc(CHISQ)): could not find function "%>%"
x2 <- fread(files[2]) %>% arrange(P) # repeat it until no more independent signal can be found.
#> Error in fread(files[2]) %>% arrange(P): could not find function "%>%"
x3 <- fread(files[3]) %>% arrange(P)
#> Error in fread(files[3]) %>% arrange(P): could not find function "%>%"
x4 <- fread(files[4]) %>% arrange(P)
#> Error in fread(files[4]) %>% arrange(P): could not find function "%>%"
x5 <- fread(files[5]) %>% arrange(P)
#> Error in fread(files[5]) %>% arrange(P): could not find function "%>%"
x6 <- fread(files[6]) %>% arrange(P)
#> Error in fread(files[6]) %>% arrange(P): could not find function "%>%"
x7 <- fread(files[7]) %>% arrange(P)
#> Error in fread(files[7]) %>% arrange(P): could not find function "%>%"
x8 <- fread(files[8]) %>% arrange(P)
#> Error in fread(files[8]) %>% arrange(P): could not find function "%>%"

env <- ls() # get the environment; this line if were put in main func will lead to error.
environment.list <- tibble(item = as.vector(grep("^x[0-9]+$", x = env, value = TRUE))) # load the environment
#> Error in tibble(item = as.vector(grep("^x[0-9]+$", x = env, value = TRUE))): could not find function "tibble"
check_significant_SNP(x8, environment.list, significance_level = 5e-8)
#> Error: object 'environment.list' not found

# You can mannually check the p-value of one SNP in previous environment.list
```
