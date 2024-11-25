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
#' @importFrom viridis scale_color_viridis
#' @import ggplot2
leo_scale_color <- function(plot, color_palette = "npg") {
  # Check if color_palette is a recognized ggsci palette
  if (color_palette == "npg") {
    plot <- plot + ggsci::scale_color_npg()
  } else if (color_palette == "lancet") {
    plot <- plot + ggsci::scale_color_lancet()
  } else if (color_palette == "jama") {
    plot <- plot + ggsci::scale_color_jama()
  } else if (color_palette == "nejm") {
    plot <- plot + ggsci::scale_color_nejm()
  } else if (color_palette == "d3") {
    plot <- plot + ggsci::scale_color_d3()
  } else if (color_palette == "tron") {
    plot <- plot + ggsci::scale_color_tron()
  } else if (color_palette == "igv") {
    plot <- plot + ggsci::scale_color_igv()
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
  # Check if color_palette is a viridis palette
  else if (color_palette %in% c("magma", "inferno", "plasma", "viridis", "cividis")) {
    plot <- plot + viridis::scale_color_viridis(option = color_palette)
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
#' This function creates a scatter plot to visualize the correlation between two vectors,
#' displaying the correlation coefficient and p-value on the plot.
#' It uses `ggplot2` for visualization and `ggsci` for color schemes.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs annotate theme_minimal theme element_text element_blank
#'
#' @param vector_x A numeric vector.
#' @param vector_y A numeric vector of the same length as \code{vector_x}.
#' @param method A character string specifying the correlation method ("spearman" or "pearson").
#'               Defaults to "spearman".
#' @param color_palette A character string or vector specifying the color palette to use.
#'               Can be a palette name from `ggsci`, `RColorBrewer`, `viridis`, or a custom color vector.
#' @param title A character string for the plot title. Defaults to "Correlation Plot".
#' @param xlab A character string for the x-axis label. Defaults to "Vector X".
#' @param ylab A character string for the y-axis label. Defaults to "Vector Y".
#' @param point_size Numeric value for point size in the scatter plot. Defaults to 4.
#' @param point_color A character string specifying color for the scatter points. Defaults to "#BB7CD8".
#' @param alpha Numeric value for the transparency of points. Defaults to 0.75.
#' @param line_color A character string specifying color for the trend line. Defaults to "#BB7CD8".
#' @param line_type Character specifying the type of line ("solid", "dashed", etc.). Defaults to "dashed".
#' @param line_size Numeric value specifying the thickness of the trend line. Defaults to 1.2.
#' @param ci_alpha Numeric value for the transparency level of the confidence interval. Defaults to 0.6.
#' @param title_size Numeric value for title text size. Defaults to 14.
#' @param xlab_size Numeric value for x-axis label text size. Defaults to 12.
#' @param ylab_size Numeric value for y-axis label text size. Defaults to 12.
#' @param axis_text_size Numeric value for axis text size. Defaults to 10.
#' @param ... Additional arguments passed to \code{\link[stats]{cor.test}}.
#' @param point_stroke It could also be NA!!!
#'
#' @return A ggplot object representing the scatter plot with correlation information.
#' @export
#' @seealso \code{\link{correlation_calculate}} for calculating correlation coefficients and p-values.
#' @seealso [leo_scale_color()] for applying color palettes to ggplot objects.
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
