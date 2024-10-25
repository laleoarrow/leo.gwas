#' Convert CHR:BP to rsID
#'
#' This function takes a dataframe that must include columns 'CHR' and 'BP',
#' and it appends the corresponding rsID by querying the SNPlocs.Hsapiens.dbSNP155.GRCh37 database.
#' The function returns a dataframe with the rsIDs included.
#'
#' @param df A dataframe containing at least the columns 'CHR' and 'BP'.
#' @param ref str indicating ref version.
#' - "GRCh37" for GRCh37 Ref panel
#' - "GRCh38" for GRCh38 Ref panel
#' @return A dataframe with an additional 'RefSNP_id' column which contains the rsIDs.
#' @importFrom data.table ":="
#' @examples
#' pacman::p_load(data.table, BSgenome, leo.gwas)
#' library("SNPlocs.Hsapiens.dbSNP155.GRCh37") # for GRCh37
#' library("SNPlocs.Hsapiens.dbSNP155.GRCh38") # for GRCh38
#' df <- data.frame(
#'   CHR = c(1, 1),
#'   BP = c(15211, 15820)
#' )
#' result <- add_rsid(df); result
#' @export
add_rsid <- function(dat, ref = "GRCh37") {
  if (!("CHR" %in% colnames(dat)) || !("BP" %in% colnames(dat))) {
    stop("DataFrame must contain 'CHR' and 'BP' columns")
  }

  # Convert dat to data.table if it is not one already
  # library(data.table)
  if (!data.table::is.data.table(dat)) {data.table::setDT(dat)}
  dat[, ranges := paste0(CHR, ":", BP, "-", BP)]

  # Load SNP data - assuming GRCh37, modify if using GRCh38
  if (ref == "GRCh37") {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37
  } else {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
  }

  # Find overlaps
  message(paste0("Translating RSID using: \n - BSgenome::snpsByOverlaps() \n - With ", snps@data_pkgname))
  snp.res <- BSgenome::snpsByOverlaps(snps, GRanges(dat$ranges))

  # Convert results to data.table
  snp.res.dt <- as.data.table(snp.res)
  snp.res.dt[, ranges := stringr::str_c(seqnames, ":", pos, "-", pos)]

  # Merge data
  trans.dat <- merge(snp.res.dt, dat, by = "ranges")
  columns_to_remove <- c("strand", "alleles_as_ambig", "ranges", "seqnames", "pos")
  trans.dat[, (columns_to_remove) := NULL]; gc()

  # Drop NA rsid and return
  # trans.dat <- trans.dat %>% drop_na(RefSNP_id)
  leo.gwas::leo_message("Remember to check if there is any NA in the RefSNP_id column.")
  leo.gwas::leo_message(">>> table(is.na(dat$RefSNP_id))")
  leo.gwas::leo_message(">>> vkh_meta %>% map_dbl(~sum(is.na(.)))")

  return(trans.dat)
}

#' Add Chromosome Position (CHR and POS) to Dataset Based on rsID
#'
#' This function takes a dataset with a column containing rsIDs (SNP IDs) and
#' adds the corresponding chromosome (CHR) and position (POS) information.
#' It queries the SNPlocs.Hsapiens.dbSNP155.GRCh37 database (or GRCh38 if specified) to retrieve the genomic positions.
#' The function returns a dataframe with the additional 'CHR' and 'POS' columns appended.
#'
#' @param dat A dataframe containing at least a column with SNP IDs (rsIDs).
#' @param snp_col A string indicating the column name containing SNP IDs (default is "SNP").
#' @param ref A string indicating the reference genome version. Default is "GRCh37", can also use "GRCh38".
#' @return A dataframe with additional 'CHR' and 'POS' columns.
#' @importFrom data.table ":="
#' @examples
#' pacman::p_load(data.table, dplyr, BSgenome, SNPlocs.Hsapiens.dbSNP155.GRCh37)
#' zuo_ref <- fread("/path/to/1KG-EAS-EAF.txt.gz") # Input dataset with rsID (SNP column)
#' result <- add_chrpos(zuo_ref, snp_col = "SNP", ref = "GRCh37")
#' @export
add_chrpos <- function(dat, snp_col = "SNP", ref = "GRCh37") {

  # Check if the SNP column exists in the dataframe
  if (!(snp_col %in% colnames(dat))) {
    stop(paste0("The column '", snp_col, "' is not found in the dataframe."))
  }

  # Load SNP reference data
  if (ref == "GRCh37") {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37
    leo.gwas::leo_message("ℹ Loading SNPlocs.Hsapiens.dbSNP155.GRCh37 database.")
  } else if (ref == "GRCh38") {
    snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
    leo.gwas::leo_message("ℹ Loading SNPlocs.Hsapiens.dbSNP155.GRCh38 database.")
  } else {
    stop("Invalid reference version. Please choose 'GRCh37' or 'GRCh38'.")
  }

  # Get CHR and POS for the SNPs from the SNPlocs database
  snp_info <- data.table::setDT(data.frame(snpsById(snps, dat[[snp_col]], ifnotfound = "drop"))) %>%
    dplyr::select(seqnames, pos, RefSNP_id) %>%
    dplyr::rename(CHR = seqnames, POS = pos, SNP = RefSNP_id)

  # Merge the SNP info with the original dataset, keeping all original columns
  merged_data <- merge(dat, snp_info, by.x = snp_col, by.y = "SNP", all.x = TRUE) %>%
    select(SNP, CHR, POS, everything())

  # Check for missing SNPs (those not found in the SNP reference)
  missing_snps <- merged_data[is.na(merged_data$CHR), ]$SNP
  if (length(missing_snps) > 0) {
    message(paste0("Cannot find the following SNPs:\n", paste(missing_snps, collapse = ", "), "\n"))
  }

  # Return the final dataset with CHR and POS columns added
  return(merged_data)
}

#' Map Gene Symbols to Genomic Positions
#'
#' This function maps gene symbols to their genomic positions (chromosome,
#' start, end, strand) using the specified method and genome assembly.
#'
#' The function supports two methods:
#' - `"bioconductor"`: uses Bioconductor packages. See [map_using_bioconductor()].
#' - `"gtf"`: uses a GTF file. See [map_using_gtf()].
#'
#' @param genes A character vector of gene symbols or a data frame containing gene symbols.
#' @param gene_col The column name of gene symbols if `genes` is a data frame.
#' @param method The method to use: `"bioconductor"` or `"gtf"`.
#' @param genome The genome assembly to use: `"hg19"` or `"hg38"`.
#' @param ... Additional arguments to pass to the GTF method.
#'
#' @return A data frame with columns: gene_symbol, chr, bp_start, bp_end, strand.
#' @importFrom dplyr %>% select mutate left_join filter pull
#' @importFrom rlang sym
#' @examples
#' \dontrun{
#' # Using Bioconductor method with a character vector of gene symbols
#' leo_map_gene_to_chrbp(genes = c("TP53", "BRCA1", "EGFR"), method = "bioconductor", genome = "hg19")
#'
#' # Using Bioconductor method with a data frame
#' leo_map_gene_to_chrbp(genes = data.frame(gene_name = c("TP53", "BRCA1", "EGFR"), value = c(1.2, 3.4, 5.6)),
#'                       gene_col = "gene_name", method = "bioconductor", genome = "hg19")
#'
#' # Using GTF method with a character vector of gene symbols
#' leo_map_gene_to_chrbp(genes = c("TP53", "BRCA1", "EGFR"), method = "gtf", genome = "hg38")
#'
#' # Using GTF method with a data frame
#' leo_map_gene_to_chrbp(genes = data.frame(gene_name = c("TP53", "BRCA1", "EGFR"), value = c(1.2, 3.4, 5.6)),
#'                       gene_col = "gene_name", method = "gtf", genome = "hg38", download_dir = "~/project/ref/gtf")
#' }
#' @export
leo_map_gene_to_chrbp <- function(genes,
                                  gene_col = NULL,
                                  method = c("bioconductor", "gtf"),
                                  genome = c("hg19", "hg38"),
                                  ...) {
  method <- match.arg(method)
  genome <- match.arg(genome)

  if (method == "bioconductor") {
    result <- map_using_bioconductor(genes = genes, gene_col = gene_col, genome = genome)
  } else if (method == "gtf") {
    result <- map_using_gtf(genes = genes, gene_col = gene_col, genome = genome, ...)
  } else {
    stop("Invalid method specified.")
  }

  return(result)
}


#' Map Gene Symbols Using Bioconductor Packages
#'
#' @param genes A character vector of gene symbols or a data frame containing gene symbols.
#' @param gene_col The column name of gene symbols if `genes` is a data frame.
#' @param genome The genome assembly to use: `"hg19"` or `"hg38"`.
#'  - For hg19, it needs `TxDb.Hsapiens.UCSC.hg19.knownGene` Bioconductor package
#'  - For hg38, it needs `TxDb.Hsapiens.UCSC.hg38.knownGene` Bioconductor package
#' @return A data frame with mapping results.
#' @importFrom dplyr %>% mutate left_join select pull
#' @importFrom rlang sym
#' @importFrom AnnotationDbi mapIds
#' @importFrom GenomicFeatures genes
#' @examples
#' \dontrun{
#' gene_symbols <- c("TP53", "BRCA1", "EGFR") # Example with gene symbol vector
#' map_using_bioconductor(genes = gene_symbols, genome = "hg19")
#'
#' gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInformation = c(1,2,3)) # Example with data frame input
#' map_using_bioconductor(genes = gene_symbols_df, gene_col = "GeneName" , genome = "hg19")
#' }
#' @export
map_using_bioconductor <- function(genes, gene_col = NULL, genome = c("hg19", "hg38")) {
  genome <- match.arg(genome)

  # Load the appropriate TxDb package based on genome assembly
  if (genome == "hg19") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
      message("Package 'TxDb.Hsapiens.UCSC.hg19.knownGene' is required. You can install it using BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene').")
      stop("Package 'TxDb.Hsapiens.UCSC.hg19.knownGene' is required.")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if (genome == "hg38") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      message("Package 'TxDb.Hsapiens.UCSC.hg38.knownGene' is required. You can install it using BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene').")
      stop("Package 'TxDb.Hsapiens.UCSC.hg38.knownGene' is required.")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }

  # Load gene annotation package
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("Package 'org.Hs.eg.db' is required. You can install it using BiocManager::install('org.Hs.eg.db').")
    stop("Package 'org.Hs.eg.db' is required.")
  }

  # Process input genes
  if (is.data.frame(genes)) {
    if (is.null(gene_col) || !(gene_col %in% colnames(genes))) {
      stop("Please specify a valid 'gene_col' that exists in the data frame.")
    }
    gene_col_sym <- rlang::sym(gene_col)
    gene_symbols <- genes %>% dplyr::pull(!!gene_col_sym)
    input_df <- genes
  } else if (is.vector(genes)) {
    gene_symbols <- genes
    input_df <- data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE)
  } else {
    stop("Input 'genes' must be a vector or a data frame.")
  }

  # Map gene symbols to Entrez IDs
  entrez_ids <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_symbols,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
  )

  # Retrieve genomic positions using Entrez IDs
  suppressMessages( gene_locations <- GenomicFeatures::genes(txdb, columns = c("gene_id"), single.strand.genes.only = T) )
  gene_locations_df <- as.data.frame(gene_locations) %>%
    dplyr::mutate(gene_id = as.character(gene_id))

  # Merge & Mapping
  mapping_df <- data.frame(
    gene_symbol = gene_symbols,
    entrez_id = entrez_ids,
    stringsAsFactors = FALSE
  )

  result_df <- mapping_df %>%
    dplyr::left_join(gene_locations_df, by = c("entrez_id" = "gene_id")) %>%
    dplyr::select(gene_symbol, chr = seqnames, bp_start = start, bp_end = end, strand) %>%
    dplyr::mutate(chr = gsub("^chr", "", chr))

  # If input was a data frame, append the results
  if (is.data.frame(genes)) {
    final_df <- genes %>% dplyr::left_join(result_df, by = stats::setNames("gene_symbol", gene_col))
    return(final_df)
  } else {
    return(result_df)
  }
}

#' Map Gene Symbols Using GTF File
#'
#' @param genes A character vector of gene symbols or a data frame containing gene symbols.
#' @param gene_col The column name of gene symbols if `genes` is a data frame.
#' @param genome The genome assembly to use: `"hg19"` or `"hg38"`.
#' @param gtf_file The path to a GTF file. If `NULL`, the function will download the appropriate GTF file.
#' @param download_dir The path where you wanna store the downloaded gtf file
#'
#' @return A data frame with mapping results.
#' @importFrom dplyr %>% mutate left_join select filter pull
#' @importFrom rlang sym
#' @importFrom GenomicFeatures makeTxDbFromGFF genes
#' @importFrom rtracklayer import
#' @examples
#' \dontrun{
#' gene_symbols <- c("TP53", "BRCA1", "EGFR")
#' map_using_gtf(genes = gene_symbols, genome = "hg38")
#'
#' gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInformation = c(1,2,3))
#' map_using_gtf(genes = gene_symbols_df, gene_col = "GeneName" , genome = "hg19")
#' }
#' @export
map_using_gtf <- function(genes, gene_col = NULL, genome = c("hg19", "hg38"), gtf_file = NULL, download_dir = "~/project/ref/gtf") {
  genome <- match.arg(genome)

  # Process input genes
  if (is.data.frame(genes)) {
    if (is.null(gene_col) || !(gene_col %in% colnames(genes))) {
      stop("Please specify a valid 'gene_col' that exists in the data frame.")
    }
    gene_col_sym <- rlang::sym(gene_col)
    gene_symbols <- genes %>% dplyr::pull(!!gene_col_sym)
    input_df <- genes
  } else if (is.vector(genes)) {
    gene_symbols <- genes
    input_df <- data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE)
  } else {
    stop("Input 'genes' must be a vector or a data frame.")
  }

  # Set default download directory
  if (is.null(download_dir)) { download_dir <- "~/project/ref/gtf" }
  if (!dir.exists(download_dir)) { dir.create(download_dir, recursive = TRUE) }
  message(paste0(">>> Setting the download_dir: ", download_dir, "\n"))

  # Download GTF file if not provided
  if (is.null(gtf_file)) {
    if (genome == "hg19") {
      # gtf_url <- "ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
      gtf_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz"
      gtf_file <- file.path(download_dir, basename(gtf_url))
    } else if (genome == "hg38") {
      # gtf_url <- "ftp://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz"
      gtf_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz"
      gtf_file <- file.path(download_dir, basename(gtf_url))
    }
    if (!file.exists(gtf_file)) {
      message("Downloading GTF file...")
      utils::download.file(gtf_url, destfile = gtf_file)
    }
  }

  # Create TxDb object from GTF file
  message("Creating TxDb object from GTF file...")
  txdb <- suppressWarnings( GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf") )

  # Retrieve gene locations
  gene_locations <- GenomicFeatures::genes(txdb)
  gene_locations_df <- as.data.frame(gene_locations)

  # Import GTF file to get gene symbols
  message("Importing GTF file to get gene symbols...")
  gtf_data <- rtracklayer::import(gtf_file)

  # Determine available attribute columns
  available_cols <- names(S4Vectors::mcols(gtf_data))
  if ("gene_name" %in% available_cols) {
    gene_id_col <- "gene_id"
    gene_name_col <- "gene_name"
  } else if ("gene" %in% available_cols) {
    gene_id_col <- "gene_id"
    gene_name_col <- "gene"
  } else {
    stop("Neither 'gene_name' nor 'gene' columns are present in the GTF file.")
  }
  gene_symbols_map <- unique(S4Vectors::mcols(gtf_data)[, c(gene_id_col, gene_name_col)])
  colnames(gene_symbols_map) <- c("gene_id", "gene_name")
  gene_symbols_map <- as.data.frame(gene_symbols_map)

  # Merge gene locations with gene symbols
  gene_locations_df <- gene_locations_df %>%
    dplyr::left_join(gene_symbols_map, by = c("gene_id"))

  # Filter for input gene symbols
  result_df <- gene_locations_df %>%
    dplyr::filter(gene_name %in% gene_symbols) %>%
    dplyr::select(gene_symbol = gene_name, chr = seqnames, bp_start = start, bp_end = end, strand)

  # Clean up chromosome names and strand
  result_df <- result_df %>%
    dplyr::mutate(
      chr = sub("^NC_0*([0-9XYMT]+)\\..*$", "\\1", chr)
    )

  # If input was a data frame, append the results
  if (is.data.frame(genes)) {
    final_df <- genes %>%
      dplyr::left_join(result_df, by = stats::setNames("gene_symbol", gene_col))
    return(final_df)
  } else {
    return(result_df)
  }
}
