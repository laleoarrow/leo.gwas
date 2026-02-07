# Apply Color Palette to ggplot

This function applies a specified color palette to a ggplot object,
supporting `ggsci`, `RColorBrewer`, and `viridis` palettes, as well as
custom colors.

## Usage

``` r
leo_scale_color(plot, color_palette = "npg")
```

## Arguments

- plot:

  A ggplot object to which the color palette will be applied.

- color_palette:

  A character string specifying the color palette to use. Options
  include `ggsci` palettes ("npg", "lancet", "jama", etc.),
  `RColorBrewer` palettes, `viridis` palettes, or a custom vector of
  colors.

## Value

A ggplot object with the applied color palette.
