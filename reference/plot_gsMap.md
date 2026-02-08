# Plot gsMap (full & highlight) with robust color mapping

Main plotting function. Full panel uses per-annotation colors; highlight
panel keeps all data but greys out non-highlighted levels.

## Usage

``` r
plot_gsMap(
  path,
  width1 = 6,
  height1 = 4,
  width2 = 8,
  height2 = 4,
  color = NULL,
  anno_colors = NULL,
  reverse_x = FALSE,
  reverse_y = TRUE,
  save_folder = "./figure/gsmap/tmp",
  basename = NULL,
  traitname = NULL,
  highlight_tissue = NULL,
  highlight_color = NULL,
  other_grey = "grey60"
)
```

## Arguments

- path:

  CSV file path with columns: sx, sy, annotation, logp.

- width1, height1:

  PDF size for full-annotation figure.

- width2, height2:

  PDF size for combined highlight+logp figure.

- color:

  Named color vector for annotations (fallback).

- anno_colors:

  Preferred named color vector for annotations (partial allowed).

- reverse_x, reverse_y:

  Reverse axes or not.

- save_folder:

  Output folder.

- basename:

  Optional plot basename.

- traitname:

  Optional trait name.

- highlight_tissue:

  Levels to highlight (default: first level if missing).

- highlight_color:

  Highlight color (single string or named vector).

- other_grey:

  Color used for non-highlight levels in highlight panel.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate & write a small CSV
set.seed(1)
df <- data.frame(sx = rnorm(60), sy = rnorm(60),
                 annotation = rep(c("Heart","Liver","Lung","Brain"), each = 15),
                 logp = runif(60, 1, 12))
fp <- file.path(tempdir(), "A_B_trait_gsMap_plot.csv")
write.csv(df, fp, row.names = FALSE)

# Partial colors; highlight Brain in gold; others grey
anno_cols <- c(Heart = "#E64B35FF", Liver = "#4DBBD5FF")
plot_gsMap(path = fp, width1 = 5, height1 = 3.5, width2 = 7, height2 = 3.5,
           anno_colors = anno_cols, reverse_x = FALSE, reverse_y = TRUE,
           save_folder = tempdir(), basename = "Demo.Base", traitname = "trait",
           highlight_tissue = "Brain", highlight_color = "#FFD700", other_grey = "grey60")
} # }
```
