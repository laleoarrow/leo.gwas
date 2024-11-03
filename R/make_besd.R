#' This contains need-to-have functions for making BESD files for SMR and HEIDI



#' Filter Chromosomes Based on SNP P-value Threshold
#'
#' This function adds an index column to the input data frame and filters chromosomes based on whether any SNP within the chromosome crosses a specified threshold. If a chromosome has at least one SNP that meets the threshold, all SNPs on that chromosome are retained. Otherwise, all SNPs on that chromosome are removed.
#'
#' @param df A data frame containing SNP data.
#' @param chr_col Character string specifying the name of the chromosome column. Default is `"CHR"`.
#' @param snp_col Character string specifying the name of the SNP identifier column. Default is `"Variant_ID"`.
#' @param p_val_col Character string specifying the name of the p-value column. Default is `"nominal_P_value"`.
#' @param threshold Numeric value specifying the p-value threshold. Default is `5e-8`.
#'
#' @return A filtered data frame with an added `index` column.
#' @export
#' @examples
#' library(dplyr)
#' eqtl_data <- data.frame(
#' CHR = c("1", "1", "2", "2", "3"),
#' Variant_ID = c("rs1", "rs2", "rs3", "rs4", "rs5"),
#' nominal_P_value = c(1e-9, 0.05, 0.2, 1e-7, 0.3)
#' );eqtl_data
#' df_filtered <- filter_chr_basedonSNP_p(
#'   df = eqtl_data,
#'   chr_col = "CHR",
#'   snp_col = "Variant_ID",
#'   p_val_col = "nominal_P_value",
#'   threshold = 5e-8
#' );df_filtered
filter_chr_basedonSNP_p <- function(df,
                               chr_col = "CHR",
                               snp_col = "Variant_ID",
                               p_val_col = "nominal_P_value",
                               threshold = 1.57e-3) {
  chr_sym <- sym(chr_col)
  snp_sym <- sym(snp_col)
  pval_sym <- sym(p_val_col)

  df_filtered <- df %>%
    dplyr::group_by(!!chr_sym) %>%
    dplyr::filter(any((!!pval_sym) < threshold)) %>%
    dplyr::ungroup()
  return(df_filtered)
}

#' Filter Chromosomes Based on SNP P-value Threshold for `.qtltoolsnomi` files
#'
#' This function adds an index column to the input data frame and filters chromosomes based on whether any SNP within the chromosome crosses a specified threshold. If a chromosome has at least one SNP that meets the threshold, all SNPs on that chromosome are retained. Otherwise, all SNPs on that chromosome are removed.
#'
#' @param df A data frame containing SNP data.
#' @param chr_col Character string specifying the name of the chromosome column. Default is `"CHR"`.
#' @param snp_col Character string specifying the name of the SNP identifier column. Default is `"Variant_ID"`.
#' @param gene_col Character string specifying the name of the gene column. Default is `"Gene"`.
#' @param p_val_col Character string specifying the name of the p-value column. Default is `"nominal_P_value"`.
#' @param threshold Numeric value specifying the p-value threshold. Default is `5e-8`.
#'
#' @return A filtered data frame with an added `index` column.
#' @export
#' @examples
#' library(dplyr)
#' eqtl_data <- data.frame(
#' Gene = c("G1", "G1", "G2", "G2", "G3"),
#' CHR = c("1", "1", "2", "2", "3"),
#' Variant_ID = c("rs1", "rs2", "rs3", "rs4", "rs5"),
#' nominal_P_value = c(1e-9, 0.05, 0.2, 1e-7, 0.3)
#' );eqtl_data
#' df_filtered <- filter_chr_basedonSNP_p_qtltools(
#'   df = eqtl_data,
#'   chr_col = "CHR",
#'   snp_col = "Variant_ID",
#'   p_val_col = "nominal_P_value",
#'   threshold = 5e-8
#' );df_filtered
filter_chr_basedonSNP_p_qtltools <- function(df,
                                    chr_col = "CHR",
                                    snp_col = "Variant_ID",
                                    gene_col = "Gene",
                                    p_val_col = "nominal_P_value",
                                    threshold = 1.57e-3) {
  chr_sym <- sym(chr_col)
  snp_sym <- sym(snp_col)
  pval_sym <- sym(p_val_col)
  gene_sym <- sym(gene_col)

  df_filtered <- df %>%
    dplyr::group_by(!!gene_sym, !!chr_sym) %>%
    dplyr::filter(any((!!pval_sym) < threshold)) %>%
    dplyr::ungroup()
  return(df_filtered)
}
