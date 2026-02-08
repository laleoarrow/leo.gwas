# Build color maps for gsMap plots (helper function)

gsMap_plot Helper to build color mappings (partial fill, highlight
override).

## Usage

``` r
plot_gsMap_color(
  annos,
  anno_colors = NULL,
  color = NULL,
  highlight_tissue = NULL,
  highlight_color = NULL
)
```

## Arguments

- annos:

  Character vector of annotation levels (unique, ordered).

- anno_colors:

  Named colors for some/all levels (preferred over `color`).

- color:

  Named colors as fallback if `anno_colors` is NULL.

- highlight_tissue:

  Character vector of levels to highlight; fallback to first of `annos`
  if none valid.

- highlight_color:

  Single color string or a named vector mapping highlight levels to
  colors.

## Value

A list with:

- final_map:

  named colors for all `annos` after auto-fill & highlight override (for
  full plot)

- hi_levels:

  chosen highlight levels

- hi_map:

  named colors for highlight levels only (for highlight-only plot)

## Examples

``` r
if (FALSE) { # \dontrun{
library(patchwork)
library(ggplot2)
# Simulate a small dataframe
df <- data.frame(sx = rnorm(30), sy = rnorm(30),
                 annotation = rep(c("TissueA","TissueB","TissueC"), each = 10),
                 logp = runif(30, 1, 10))

# User provides partial colors; others will be auto-filled
cm <- plot_gsMap_color(annos = unique(df$annotation),
                       anno_colors = c(TissueA = "#F2B701"),
                       highlight_tissue = "TissueB",
                       highlight_color  = "#FF0000")

# Full plot (manual scale with final_map)
gg <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = annotation)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = cm$final_map)

# Highlight-only plot (only highlight levels mapped)
gg_hi <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = annotation)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = cm$hi_map, name = "Highlight")
gg | gg_hi
} # }
```
