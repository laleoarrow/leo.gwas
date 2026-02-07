# Combine `.fdr` files for one or multiple outcomes (seperately)

Combine `.fdr` files for one or multiple outcomes (seperately)

## Usage

``` r
combine_smr_res_1outcome(
  dir,
  outcome,
  out_dir = file.path(dir, "combine_1outcome")
)
```

## Arguments

- dir:

  Character. The parent folder that contains subfolders with fdr files.

- outcome:

  Character or character vector. One or more outcomes to search in file
  names.

- out_dir:

  Character. Output folder; default is to create "combine_1outcome"
  under `dir`.

## Value

         NULL. This function writes combined files to disk directly.

## Examples

``` r
if (FALSE) { # \dontrun{
combine_smr_res_1outcome("/Users/leoarrow/project/iridocyclitis/output/smr", "iri3")
} # }
```
