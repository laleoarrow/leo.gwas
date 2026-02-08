##################################################
# This script is for functions for visualization #
# @author: Lu Ao (luao@stu.cqmu.edu.cn)          #
##################################################
#  ------------------------------ Basics Function ------------------------------
#' Apply Color Palette to ggplot
#'
#' This function applies a specified color palette to a ggplot object,
#' supporting `ggsci`, `RColorBrewer`, and `viridis` palettes, as well as custom colors.
#'
#' @param plot A ggplot object to which the color palette will be applied.
#' @param color_palette A character string specifying the color palette to use.
#'                      Options include `ggsci` palettes ("npg", "lancet", "jama", etc.),
#'                      `RColorBrewer` palettes, `viridis` palettes, or a custom vector of colors.
#'
#' @return A ggplot object with the applied color palette.
#' @importFrom ggsci scale_color_npg scale_color_lancet scale_color_jama scale_color_nejm scale_color_d3 scale_color_tron scale_color_igv scale_color_ucscgb scale_color_aaas scale_color_futurama scale_color_rickandmorty scale_color_simpsons
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
leo_scale_color <- function(plot, color_palette = "npg") {
  # Check if color_palette is a recognized ggsci palette
  if (color_palette == "npg") {
    plot <- plot + ggsci::scale_color_npg()
  } else if (color_palette == "lancet") {
    plot <- plot + ggsci::scale_color_lancet()
  } else if (color_palette == "ucscgb") {
    plot <- plot + ggsci::scale_color_ucscgb()
  } else if (color_palette == "aaas") {
    plot <- plot + ggsci::scale_color_aaas()
  } else if (color_palette == "futurama") {
    plot <- plot + ggsci::scale_color_futurama()
  } else if (color_palette == "rickandmorty") {
    plot <- plot + ggsci::scale_color_rickandmorty()
  } else if (color_palette == "simpsons") {
    plot <- plot + ggsci::scale_color_simpsons()
  }
  # Check if color_palette is a recognized RColorBrewer palette
  else if (color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    plot <- plot + scale_color_brewer(palette = color_palette)
  }
  # Use a custom color vector if provided
  else if (is.vector(color_palette) && all(is.character(color_palette))) {
    plot <- plot + scale_color_manual(values = color_palette)
  } else {
    warning("Color palette not recognized. Default ggplot2 colors will be used.")
  }
  return(plot)
}

#' Draw Correlation between Two Vectors
#'
#' This function creates a scatter plot to visualize the correlation between two vectors, displaying the correlation coefficient and p-value on the plot.
#'
#' @param vector_x Numeric vector.
#' @param vector_y Numeric vector of the same length as \code{vector_x}.
#' @param method Correlation method: \code{"spearman"} or \code{"pearson"}. Default \code{"spearman"}.
#' @param color_palette Character scalar or vector for color palette (kept for compatibility).
#' @param title Plot title. Default \code{"Correlation Plot"}.
#' @param xlab,ylab Axis labels. Defaults \code{"Vector X"}, \code{"Vector Y"}.
#' @param point_size Point size. Default \code{1.5}.
#' @param point_color Fill color for points (shape 21). Default \code{"#BB7CD8"}.
#' @param point_stroke Numeric stroke width for point outline. If \code{NA}, treated as \code{0}. Default \code{1}.
#' @param alpha Point transparency. Default \code{0.75}.
#' @param line_color, line_type, line_size Trend line color, type, size. Defaults \code{"#BB7CD8"}, \code{"dashed"}, \code{1.2}.
#' @param ci_alpha Confidence ribbon alpha. Default \code{0.2}.
#' @param title_size,xlab_size,ylab_size,axis_text_size Font sizes. Defaults \code{16}, \code{14}, \code{14}, \code{14}.
#' @param ... Additional arguments passed to \code{correlation_calculate()}.
#' @return A \pkg{ggplot2} object.
#' @export
#'
#' @seealso \code{\link{correlation_calculate}}, \code{\link{leo_scale_color}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs annotate theme_classic
#' @importFrom ggplot2 theme element_text element_blank element_line coord_cartesian
#' @importFrom rlang .data
#' @examples
#' vector_x <- c(10, 2, 3, 4, 5)
#' vector_y <- c(5, 6, 7, 8, 7)
#' correlation_draw(vector_x, vector_y, method = "pearson", point_size = 10, color_palette = "npg")
correlation_draw <- function(vector_x, vector_y, method = "spearman", color_palette = "npg",
                             title = "Correlation Plot", xlab = "Vector X", ylab = "Vector Y",
                             point_size = 1.5, point_color = "#BB7CD8", point_stroke = 1, alpha = 0.75,
                             line_color = "#BB7CD8", line_type = "dashed", line_size = 1.2,
                             ci_alpha = 0.2, title_size = 16, xlab_size = 14,
                             ylab_size = 14, axis_text_size = 14, ...) {

  # Calculate correlation
  correlation_result <- correlation_calculate(vector_x, vector_y, method = method, ...)
  correlation_coefficient <- round(correlation_result$correlation_coefficient, 3)
  p_value <- round(correlation_result$p_value, 3)

  # Create scatter plot with correlation coefficient and p-value annotation
  plot <- ggplot2::ggplot(data = data.frame(vector_x, vector_y), ggplot2::aes(x = vector_x, y = vector_y)) +
    ggplot2::geom_point(size = point_size, color = "black", fill = point_color, alpha = alpha, stroke = point_stroke, shape = 21) +
    ggplot2::geom_smooth(method = "lm", color = line_color, linetype = line_type,
                         size = line_size, se = TRUE, fill = line_color, alpha = ci_alpha) +  # Add confidence interval shading
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = paste("Correlation:", correlation_coefficient, "\nP-value:", p_value),
      hjust = 1.1, vjust = 1.2,
      size = 5, color = "black"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size),
      axis.text = ggplot2::element_text(size = axis_text_size, color = "black"),
      axis.title.x = ggplot2::element_text(size = xlab_size, color = "black"),
      axis.title.y = ggplot2::element_text(size = ylab_size, color = "black"),
      axis.line = ggplot2::element_line(color = "black"),      # Set axis line color
      axis.ticks = ggplot2::element_line(color = "black"),     # Set axis ticks color
      panel.grid = ggplot2::element_blank()                    # Remove grid
    )

  # Apply color palette using the helper function
  # plot <- leo_scale_color(plot, color_palette)

  return(plot)
}

#' Draw Group Comparison for a Continuous Variable
#'
#' This function creates a violin + boxplot to visualize the distribution of a continuous variable
#' across groups of a categorical variable, with optional jitter points and annotated p-value.
#'
#' @param df A data frame (optional if vectors provided).
#' @param x_col Name of categorical variable in df.
#' @param y_col Name of numeric variable in df.
#' @param vector_x Optional categorical vector.
#' @param vector_y Optional numeric vector.
#' @param test_method Statistical test: \code{"wilcox"} (default, Mann–Whitney U test) or \code{"t.test"}.
#' @param title Plot title.
#' @param alpha Transparency for violin, boxplot and jitter.
#' @param violin_fill Fill color for violin.
#' @param box_fill Fill color for boxplot.
#' @param box_color Outline color for boxplot.
#' @param jitter Logical, whether to show jitter. Default TRUE.
#' @param jitter_size Size of jitter points. Default \code{2}.
#' @param jitter_color Color of jitter points. Default \code{"black"}.
#' @param xlab,ylab Axis labels.
#' @param title_size Font size for plot title. Default \code{16}.
#' @param xlab_size Font size for x-axis label. Default \code{14}.
#' @param ylab_size Font size for y-axis label. Default \code{14}.
#' @param axis_text_size Font size for axis text. Default \code{14}.
#' @param drop_na Logical, whether to drop NA group from x variable. Default \code{FALSE}.
#' @param annotate_n Logical, whether to annotate sample size N for each group. Default \code{TRUE}.
#' @param main_color Easy way to set the vibe. Good luck trying!
#'
#' @return A ggplot2 object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_jitter labs annotate theme_classic theme element_text
#' @importFrom dplyr %>%
#' @examples
#' df <- tibble::tibble(group = rep(c(0,1), each=50), prs = rnorm(100))
#' group_comparison_draw(df, "group", "prs")
#' group_comparison_draw(df, "group", "prs", annotate_n = TRUE)
#' group_comparison_draw(vector_x = rep(c(0,1), each=50), vector_y = rnorm(100))
#' group_comparison_draw(vector_x = rep(c(0,1), each=50), vector_y = rnorm(100), jitter = FALSE)
group_comparison_draw <- function(df, x_col, y_col, vector_x = NULL, vector_y = NULL,
                                  test_method = "wilcox", title = "Group Comparison",
                                  xlab = NULL, ylab = NULL, alpha = 0.8,
                                  main_color = "#BB7CD8",
                                  violin_fill = main_color,
                                  box_fill = main_color, box_color = "black",
                                  jitter = T, jitter_size = 2, jitter_color = "black",
                                  drop_na = F, annotate_n = T,
                                  title_size = 16, xlab_size = 14, ylab_size = 14, axis_text_size = 14) {
  # prepare data
  if (!is.null(vector_x) & !is.null(vector_y)) {
    df <- data.frame(x = vector_x, y = vector_y)
    x_col <- "x"; y_col <- "y"
  }
  if (is.null(df) | is.null(x_col) | is.null(y_col)) stop("Please provide df+x_col+y_col or vector_x+vector_y")
  if (drop_na) df <- df %>% dplyr::filter(!is.na(.data[[x_col]]))

  x <- df[[x_col]]; y <- df[[y_col]]
  if (!is.numeric(y)) stop("y_col must be numeric")

  # run test
  test_res <- if (test_method == "t.test") t.test(y ~ x) else wilcox.test(y ~ x)
  p_value <- signif(test_res$p.value, 3)

  # plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data[[x_col]]), y = .data[[y_col]])) +
    ggplot2::geom_violin(fill = violin_fill, alpha = alpha / 3, trim = F) +
    ggplot2::geom_boxplot(width = 0.2, fill = box_fill, alpha = alpha / 2, color = box_color, outlier.shape = NA)
  if (jitter) p <- p + ggplot2::geom_jitter(width = 0.15, alpha = alpha, size = jitter_size, color = jitter_color)
  if (annotate_n) {
    n_df <- df %>% dplyr::group_by(.data[[x_col]]) %>% dplyr::summarise(N = dplyr::n())
    group_levels <- levels(factor(df[[x_col]]))
    n_labels <- paste0(group_levels, "\n(N=", n_df$N[match(group_levels, n_df[[x_col]])], ")")
    p <- p + ggplot2::scale_x_discrete(labels = n_labels)
  }
  p <- p +
    ggplot2::labs(title = title,
                  x = ifelse(is.null(xlab), x_col, xlab),
                  y = ifelse(is.null(ylab), y_col, ylab)) +
    ggplot2::annotate("text", x = 1.5, y = max(y, na.rm = T),
                      label = paste0("P = ", p_value), vjust = -0.5, size = 5) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size),
                   axis.text = ggplot2::element_text(size = axis_text_size, color = "black"),
                   axis.title.x = ggplot2::element_text(size = xlab_size, color = "black"),
                   axis.title.y = ggplot2::element_text(size = ylab_size, color = "black"))
  return(p)
}

#  ------------------------------ Reginal Plot Function ------------------------------

#' Loci_plot: Calculate the LD-matrix (LD r2) for the index SNP
#' @param gwas gwas summary data that needs to select loci and calculate r2
#' @param index index snp
#' @param win window size to locally calculate the r2; set it larger than that you want to plot; # default calculate 1MB
#' @param pop only applicable under 500 snps; dont use it anyway
#' @param ld_calculation if calculate the LD locally, defaut T; use F if the index snp is a rare variant (MAF<0.01)
#' @param bfile bfile
#' @param plink_bin plinkbinr::get_plink_exe()
#'
#' @return loci data with calculated r2
#' @export
ld_ps_index <- function(gwas, index = "rs999", # ps for pre-select
                        win = 1000,
                        ld_calculation = T, bfile = "/Users/leoarrow/project/ref/1kg.v3/EAS", plink_bin = plinkbinr::get_plink_exe()){
  # extract the index SNP
  chr <- gwas %>% dplyr::filter(SNP == index) %>% pull(CHR); message(paste0("CHR for the index snp is "), chr)
  pos <- gwas %>% dplyr::filter(SNP == index) %>% pull(POS); message(paste0("POS for the index snp is "), pos)
  # extract the SNPs in the LD window
  loci <- gwas %>% dplyr::filter(CHR == chr, POS>=pos-win/2*1000, POS<=pos+win/2*1000) %>% tidyr::drop_na() %>% distinct(SNP, .keep_all = T)
  message(paste0("SNP number in the loci is "), nrow(loci))
  # calculate the LD matrix
  if (ld_calculation) {
    ld <- ieugwasr::ld_matrix(
      variants = unique(loci$SNP),
      with_alleles = F,
      # pop = pop,
      bfile = bfile,
      plink_bin = plink_bin
    )
    ld_df <- ld %>% as.data.frame()
    ld_df <- tibble::rownames_to_column(ld_df, var = "SNP")
    ld_df_index <- ld_df %>% dplyr::select(SNP, index) %>% set_names("SNP", "r") %>% mutate(r2 = r^2)
    # merge the LD matrix with the GWAS data
    loci <- loci %>% left_join(ld_df_index %>% dplyr::select(-r), by = "SNP") %>% as.data.frame()
  }
  return(loci)
}

#' Loci_plot: prepare the locus data for locuszoomr
#' @param loci_data output from ld_ps_index
#' @param gene gene loci
#' @param index_snp indexed snp
#' @param online_ld whether to use online LD; default is F
#' @param flank flank size for the locus plot
#' @return prepared data which could be pass to save_regional_plot
#' @export
locuszoomr_loc <- function(loci_data, gene, online_ld = F, index_snp, flank) {
  #   ----- loc_plot using `locuszoomr`
  loc <- locus(data = loci_data,
               gene = gene,
               index_snp = index_snp,
               chrom = "CHR", #  detect_cols(., chrom, pos, p, labs, yvar)
               pos = "POS",
               p = "P",
               labs = "SNP",
               flank = flank, # up/low flank setting; 1e5 = 100kb
               LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75")


  #  ----- LD
  if (online_ld) {loc <- link_LD(loc, pop = "EUR", token = "8866c6877cb8", method = "matrix")}
  #  ----- recomb
  library(rtracklayer)
  local_recomb_19 <- "/Users/leoarrow/project/ref/recombMap/hapMapRelease24CombinedRecombMap.bw"
  recomb.hg19 <- import.bw(local_recomb_19)
  loc <- link_recomb(loc, recomb = recomb.hg19)
  return(loc)
}

#' Loci_plot: save_regional_plot
#'
#' @param path path to store the plot; make sure the path is exist
#' @param loc output from locuszoomr_loc
#' @param gene gene
#' @param width width
#' @param save if T, will save plot to path; if F, return the plot only
#' @param labels labels; in case you need to indicate the index SNP and other SNP
#' @param border border for gene track
#' @param height height
#'
#' @export
save_regional_plot <- function(path, loc, gene, save = T, title = expression(paste(italic("CLPSL1"), " (T1D)")),
                               labels = c("index"), filter_gene_biotype = c("protein_coding"), border = F, width = 7.5, height = 5.5){
  # Check if the path exists; interactively create the path if not
  if (!dir.exists(dirname(path))) {
    create_dir <- readline(prompt = "Directory does not exist. Do you want to create it? (yes/no): ")
    if (tolower(create_dir) == "yes") { # i.e., input case in-sensitive
      dir.create(dirname(path), recursive = TRUE)
      message(paste("Directory created >>>", dirname(path)))
    } else {
      stop("Directory does not exist and was not created.")
    }
  }
  # Check if the plot already exists to prevent overwriting it because you forget to change the path
  if (file.exists(path)) {
    file.exist.status <- readline(prompt = "Plot already exists! Do you mean to overwrite it or forget to change the path? (yes/no): ")
    if (tolower(file.exist.status) == "yes") { # i.e., input case in-sensitive
      message("Ok, processing the plot...")
    } else {
      stop("Process aborted now.")
    }
  }
  # Save the plot
  if (save) {message(paste("Plot save to >>>", path));pdf(path, width = width, height = height)}
  # plot
  locus_plot(loc, legend_pos = "topright",
             labels = labels,
             filter_gene_biotype = filter_gene_biotype,
             highlight_col = "#E64B35FF",
             use_layout = T,
             border = border,
             recomb_col = "#4DBBD5FF",
             highlight = gene)
  title(main=title)
  if (save) {dev.off()}
}

# ------------------------------ gsMap Plot Function ------------------------------
#' Plot gsMap (full & highlight) with robust color mapping
#'
#' Main plotting function. Full panel uses per-annotation colors; highlight panel keeps all data but greys out non-highlighted levels.
#'
#' @param path CSV file path with columns: sx, sy, annotation, logp.
#' @param width1,height1 PDF size for full-annotation figure.
#' @param width2,height2 PDF size for combined highlight+logp figure.
#' @param color Named color vector for annotations (fallback).
#' @param anno_colors Preferred named color vector for annotations (partial allowed).
#' @param reverse_x,reverse_y Reverse axes or not.
#' @param save_folder Output folder.
#' @param basename Optional plot basename.
#' @param traitname Optional trait name.
#' @param highlight_tissue Levels to highlight (default: first level if missing).
#' @param highlight_color Highlight color (single string or named vector).
#' @param other_grey Color used for non-highlight levels in highlight panel.
#'
#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
#' @importFrom stringr str_extract str_replace
#' @importFrom ggplot2 ggplot aes labs theme_void theme element_rect element_text guides guide_legend guide_colorbar
#' @importFrom ggplot2 scale_color_manual scale_x_reverse scale_y_reverse ggsave
#' @importFrom ggrastr geom_point_rast
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid unit
#' @importFrom glue glue
#' @importFrom leo.basic leo_log
#'
#' @examples
#' # Simulate & write a small CSV
#' set.seed(1)
#' df <- data.frame(sx = rnorm(60), sy = rnorm(60),
#'                  annotation = rep(c("Heart","Liver","Lung","Brain"), each = 15),
#'                  logp = runif(60, 1, 12))
#' fp <- file.path(tempdir(), "A_B_trait_gsMap_plot.csv")
#' write.csv(df, fp, row.names = FALSE)
#'
#' # Partial colors; highlight Brain in gold; others grey
#' anno_cols <- c(Heart = "#E64B35FF", Liver = "#4DBBD5FF")
#' plot_gsMap(path = fp, width1 = 5, height1 = 3.5, width2 = 7, height2 = 3.5,
#'            anno_colors = anno_cols, reverse_x = FALSE, reverse_y = TRUE,
#'            save_folder = tempdir(), basename = "Demo.Base", traitname = "trait",
#'            highlight_tissue = "Brain", highlight_color = "#FFD700", other_grey = "grey60")
plot_gsMap <- function(path, width1 = 6, height1 = 4, width2 = 8, height2 = 4,
                       color = NULL, anno_colors = NULL,
                       reverse_x = FALSE, reverse_y = TRUE, save_folder = "./figure/gsmap/tmp",
                       basename = NULL, traitname = NULL,
                       highlight_tissue = NULL, highlight_color = NULL,
                       other_grey = "grey60") {
  # IO & meta
  leo.basic::leo_log("Read: {path}")
  df <- vroom::vroom(path, delim = ",")
  base <- tools::file_path_sans_ext(basename(path))
  if (is.null(basename))  basename  <- stringr::str_extract(base, "^[^.]+\\.[^.]+")
  if (is.null(traitname)) traitname <- stringr::str_replace(base, "^[^_]*_[^_]*_(.*)_gsMap_plot$", "\\1")
  if (!dir.exists(save_folder)) { dir.create(save_folder, recursive = TRUE); leo.basic::leo_log("Create dir: {save_folder}") }

  # Color maps via helper
  annos <- unique(df$annotation)
  cm <- plot_gsMap_color(annos = annos, anno_colors = anno_colors, color = color,
                         highlight_tissue = highlight_tissue, highlight_color = highlight_color)

  # Full annotation plot
  leo.basic::leo_log("Plot full annotations ...")
  p0 <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = annotation)) +
    ggrastr::geom_point_rast(size = 0.1, shape = 16, raster.dpi = 400) +
    (if (reverse_x) ggplot2::scale_x_reverse() else NULL) +
    (if (reverse_y) ggplot2::scale_y_reverse() else NULL) +
    ggplot2::scale_color_manual(values = cm$final_map, name = "Annotation") +
    ggplot2::labs(title = basename, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background   = ggplot2::element_rect(fill = "black", color = NA),
                   panel.background  = ggplot2::element_rect(fill = "black", color = NA),
                   legend.background = ggplot2::element_rect(fill = "black", color = NA),
                   plot.title        = ggplot2::element_text(color = "white", size = 10, face = "bold", hjust = 0.5),
                   legend.title      = ggplot2::element_text(color = "white", size = 8),
                   legend.text       = ggplot2::element_text(color = "white", size = 8)) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 1.5)))
  out1 <- file.path(save_folder, glue::glue("{basename}_{traitname}_all_annotation.pdf"))
  ggplot2::ggsave(out1, p0, width = width1, height = height1); leo.basic::leo_log("Save: {.path {out1}}")

  # Highlight panel: keep all data; non-highlight -> grey
  leo.basic::leo_log("Plot highlight ...")
  other_map <- setNames(rep(other_grey, length(annos)), annos)
  hi_map_grey <- c(other_map, cm$hi_map); hi_map_grey <- hi_map_grey[annos]; hi_map_grey[names(cm$hi_map)] <- cm$hi_map
  p_highlight <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = annotation)) +
    ggrastr::geom_point_rast(size = 0.1, shape = 16, raster.dpi = 400) +
    (if (reverse_x) ggplot2::scale_x_reverse() else NULL) +
    (if (reverse_y) ggplot2::scale_y_reverse() else NULL) +
    ggplot2::scale_color_manual(values = hi_map_grey, name = "Highlight", breaks = cm$hi_levels) +
    ggplot2::labs(title = basename, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background     = ggplot2::element_rect(fill = "black", color = NA),
                   panel.background    = ggplot2::element_rect(fill = "black", color = NA),
                   legend.background   = ggplot2::element_blank(),
                   plot.title          = ggplot2::element_text(color = "white", size = 10, face = "bold", hjust = 0.5),
                   legend.title        = ggplot2::element_text(color = "white", size = 8, face = "bold",  hjust = 0.1),
                   legend.text         = ggplot2::element_text(color = "white", size = 8),
                   legend.position     = c(0, 0),
                   legend.justification = c(0, 0),
                   # legend.direction    = "horizontal",
                   legend.key.height = grid::unit(.6, "lines"),
                   # legend.spacing.y  = grid::unit(0.01, "cm"),
                   legend.box.margin   = ggplot2::margin(0, 0, 0, 0)) +
    ggplot2::guides(color = ggplot2::guide_legend(title.position = "top",
                                                  ncol = 1, byrow = TRUE,
                                                  override.aes = list(size = 2)))

  # logp panel
  p_logp <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = logp)) +
    ggrastr::geom_point_rast(size = 0.1, shape = 16, raster.dpi = 400) +
    ggplot2::scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                   limits = range(df$logp, na.rm = TRUE),
                                   breaks = scales::breaks_pretty(3),
                                   labels = scales::label_number(accuracy = 1),
                                   name = expression(atop(-log[10], (italic(P)~value)))) +
    (if (reverse_x) ggplot2::scale_x_reverse() else NULL) +
    (if (reverse_y) ggplot2::scale_y_reverse() else NULL) +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(plot.background   = ggplot2::element_rect(fill = "black", color = NA),
                   panel.background  = ggplot2::element_rect(fill = "black", color = NA),
                   legend.background = ggplot2::element_blank(),
                   legend.title      = ggplot2::element_text(color = "white", size = 8),
                   legend.text       = ggplot2::element_text(color = "white", size = 6),
                   legend.justification = c(0, 0),
                   legend.position   = c(0, 0),
                   legend.direction  = "horizontal") +
    ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                    barwidth  = grid::unit(1.5, "cm"),
                                                    barheight = grid::unit(0.25, "cm"),
                                                    ticks = FALSE))
  combined <- p_highlight | p_logp
  out2 <- file.path(save_folder, glue::glue("{basename}_{traitname}_highlight_logp.pdf"))
  ggplot2::ggsave(out2, combined, width = width2, height = height2); leo.basic::leo_log("Save: {.path {out2}}")
  leo.basic::leo_log("Done \u2705")
}

#' Build color maps for gsMap plots (helper function)
#'
#' gsMap_plot Helper to build color mappings (partial fill, highlight override).
#'
#' @param annos Character vector of annotation levels (unique, ordered).
#' @param anno_colors Named colors for some/all levels (preferred over `color`).
#' @param color Named colors as fallback if `anno_colors` is NULL.
#' @param highlight_tissue Character vector of levels to highlight; fallback to first of `annos` if none valid.
#' @param highlight_color Single color string or a named vector mapping highlight levels to colors.
#'
#' @return A list with:
#'   \item{final_map}{named colors for all `annos` after auto-fill & highlight override (for full plot)}
#'   \item{hi_levels}{chosen highlight levels}
#'   \item{hi_map}{named colors for highlight levels only (for highlight-only plot)}
#'
#' @importFrom leo.basic leo_log leo_discrete_color
#' @examples
#' library(patchwork)
#' library(ggplot2)
#' # Simulate a small dataframe
#' df <- data.frame(sx = rnorm(30), sy = rnorm(30),
#'                  annotation = rep(c("TissueA","TissueB","TissueC"), each = 10),
#'                  logp = runif(30, 1, 10))
#'
#' # User provides partial colors; others will be auto-filled
#' cm <- plot_gsMap_color(annos = unique(df$annotation),
#'                        anno_colors = c(TissueA = "#F2B701"),
#'                        highlight_tissue = "TissueB",
#'                        highlight_color  = "#FF0000")
#'
#' # Full plot (manual scale with final_map)
#' gg <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = annotation)) +
#'   ggplot2::geom_point() +
#'   ggplot2::scale_color_manual(values = cm$final_map)
#'
#' # Highlight-only plot (only highlight levels mapped)
#' gg_hi <- ggplot2::ggplot(df, ggplot2::aes(sx, sy, color = annotation)) +
#'   ggplot2::geom_point() +
#'   ggplot2::scale_color_manual(values = cm$hi_map, name = "Highlight")
#' gg | gg_hi
plot_gsMap_color <- function(annos, anno_colors = NULL, color = NULL,
                             highlight_tissue = NULL, highlight_color = NULL) {
  # Validate & normalize levels
  stopifnot(!is.null(annos), length(annos) > 0)
  annos <- unique(as.character(annos))

  # Choose user map (anno_colors preferred)
  user_map <- if (!is.null(anno_colors)) anno_colors else color

  # Build final_map with partial fill
  if (is.null(user_map)) {
    final_map <- leo.basic::leo_discrete_color(levels = annos); leo.basic::leo_log("Auto-colored {length(annos)} levels")
  } else {
    if (is.null(names(user_map)) || any(is.na(names(user_map)) | names(user_map) == "")) stop("anno_colors/color must be a named vector whose names match annotation levels")
    user_map <- user_map[!is.na(names(user_map)) & names(user_map) != ""]
    extra_keys <- setdiff(names(user_map), annos)
    if (length(extra_keys) > 0) {
      leo.basic::leo_log("Removed {length(extra_keys)} unknown keys: {paste(extra_keys, collapse=', ')}", level = "warning")
      user_map <- user_map[setdiff(names(user_map), extra_keys)]
    }
    missing_keys <- setdiff(annos, names(user_map))
    if (length(missing_keys) > 0) {
      auto_map <- leo.basic::leo_discrete_color(levels = missing_keys)
      final_map <- c(user_map, auto_map); leo.basic::leo_log("Filled colors for: {paste(missing_keys, collapse=', ')}")
    } else {
      final_map <- user_map; leo.basic::leo_log("Using user-provided colors")
    }
    final_map <- final_map[annos]
  }

  # Choose highlight levels
  if (!is.null(highlight_tissue) && length(highlight_tissue) > 0) {
    hi_levels <- intersect(as.character(highlight_tissue), annos)
    if (length(hi_levels) == 0) { hi_levels <- annos[1]; leo.basic::leo_log("Fallback highlight: {hi_levels}") }
  } else {
    hi_levels <- annos[1]; leo.basic::leo_log("Default highlight: {hi_levels}")
  }

  # ---- override highlight colors into final_map (robust for >1 highlights)
  if (!is.null(highlight_color)) {
    if (is.null(names(highlight_color))) {
      # No names: support two cases
      if (length(hi_levels) == 1) {
        final_map[hi_levels] <- highlight_color[1]
        leo.basic::leo_log("Highlight override: {hi_levels} -> {highlight_color[1]}")
      } else if (length(highlight_color) == length(hi_levels)) {
        # Map in order when lengths match
        final_map[hi_levels] <- highlight_color
        leo.basic::leo_log("Highlight override (ordered): {paste(hi_levels, collapse=', ')}")
      } else {
        # Length mismatch: keep original per-level colors
        leo.basic::leo_log("Multiple highlights but single/mismatched highlight_color; keep original colors. Supply a named vector or same-length vector to override.", level = "warning")
      }
    } else {
      # Named vector: override by matching names
      cover_keys <- intersect(names(highlight_color), hi_levels)
      if (length(cover_keys) == 0) {
        leo.basic::leo_log("Named highlight_color not matched; no override", level = "warning")
      } else {
        final_map[cover_keys] <- highlight_color[cover_keys]
        leo.basic::leo_log("Highlight override: {paste(paste0(cover_keys,'=',highlight_color[cover_keys]), collapse=', ')}")
      }
    }
  }

  # Build hi_map as a named vector containing only highlight levels
  hi_map <- final_map[hi_levels]
  list(final_map = final_map, hi_levels = hi_levels, hi_map = hi_map)
}
