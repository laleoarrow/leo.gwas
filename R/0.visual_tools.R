##################################################
# This script is for functions for visualization #
# @author: Lu Ao (luao@stu.cqmu.edu.cn)          #
##################################################

#  ------------------------------ Reginal Plot Function ------------------------------ #

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
#' @return
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
