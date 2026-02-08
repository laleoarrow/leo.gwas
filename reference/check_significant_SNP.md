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

## Examples

``` r
if (FALSE) { # \dontrun{
# specify the directory to store the HLA original data and subsequent
con_dir <- "/Users/leoarrow/project/VKH2024/data/zuo/con_su"
files <- list.files(con_dir,full.names = T) %>% as.vector(); files # update it each time
# and sort it by CHISQ/P value
x1 <- fread(files[1]) %>% arrange(desc(CHISQ)); head(x1) # for the 1st, just read the data 
x2 <- fread(files[2]) %>% arrange(P) # repeat until no more independent signal
x3 <- fread(files[3]) %>% arrange(P)
x4 <- fread(files[4]) %>% arrange(P)
x5 <- fread(files[5]) %>% arrange(P)
x6 <- fread(files[6]) %>% arrange(P)
x7 <- fread(files[7]) %>% arrange(P)
x8 <- fread(files[8]) %>% arrange(P)

env <- ls() # get the environment; this line if were put in main func will lead to error.
# load the environment
environment.list <- tibble(item = as.vector(grep("^x[0-9]+$", x = env, value = TRUE))) 
check_significant_SNP(x8, environment.list, significance_level = 5e-8)

# You can mannually check the p-value of one SNP in previous environment.list
} # }
```
