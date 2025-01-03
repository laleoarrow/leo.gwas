#' Convert CHR:BP to rsID (not recommended)
#'
#' This function takes a dataframe that must include columns 'CHR' and 'BP',
#' and it appends the corresponding rsID by querying the SNPlocs.Hsapiens.dbSNP155.GRCh37 database.
#' The function returns a dataframe with the rsIDs included.
#'
#' @param ref str indicating ref version.
#' - "GRCh37" for GRCh37 Ref panel
#' @param dat
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

#' Add rsID based on local reference file (recommended)
#'
#' This function takes a data frame with at least chromosome, position, and allele columns,
#' and matches rsID based on a local reference file, considering reversed/complement alleles.
#' It allows for flexibility in the names of the allele columns (e.g., 'REF'/'ALT' or 'A1'/'A2').
#' The function returns the original data frame with rsID added in the first column.
#'
#' @param dat A data frame containing at least chromosome, position, and allele columns.
#' @param local_ref A data frame containing at least 'ID' and 'SNP' columns.
#' @param chr_col Name of the chromosome column in 'dat'. Default is 'CHR'.
#' @param pos_col Name of the position column in 'dat'. Default is 'BP'.
#' @param a1_col Name of the first allele column in 'dat'.  Default is 'A1' (e.g., 'A1' or 'ALT').
#' @param a2_col Name of the second allele column in 'dat'.  Default is 'A2' (e.g., 'A2' or 'REF').
#' @return The original data frame with an added 'rsID' column as the first column.
#' @examples
#' \dontrun{
#' library(dplyr)
#' # Example data with REF and ALT
#' dat <- data.frame(
#'   CHR = c(7, 12, 4),
#'   BP = c(6013153, 126890980, 40088896),
#'   REF = c("G", "A", "T"),
#'   ALT = c("A", "G", "A")
#' )
#' # Reference data
#' local_ref <- data.frame(
#'   ID = c("7:6013153:A:G", "12:126890980:G:A", "4:40088896:A:T"),
#'   SNP = c("rs10000", "rs1000000", "rs10000000")
#' )
#' result <- add_rsid_using_ref(dat, local_ref, a1_col = "ALT", a2_col = "REF")
#' print(result)
#'
#' # Example data with A1 and A2
#' dat2 <- data.frame(
#'   CHR = c(7, 12, 4),
#'   POS = c(6013153, 126890980, 40088896),
#'   A1 = c("A", "G", "A"),
#'   A2 = c("G", "A", "T")
#' )
#' result2 <- add_rsid_using_ref(dat2, local_ref, pos_col = "POS")
#' print(result2)
#' }
#' @export
add_rsid_using_ref <- function(dat, local_ref, chr_col = "CHR", pos_col = "BP", a1_col = "A1", a2_col = "A2") {

  # Ensure required columns are present
  required_cols_dat <- c(chr_col, pos_col, a1_col, a2_col)
  required_cols_ref <- c("ID", "SNP")
  if (!all(required_cols_dat %in% colnames(dat))) {
    stop("Data frame 'dat' must contain columns: ", paste(required_cols_dat, collapse = ", "))
  }
  if (!all(required_cols_ref %in% colnames(local_ref))) {
    stop("Reference data 'local_ref' must contain columns: ", paste(required_cols_ref, collapse = ", "))
  }

  # Define function to get complement
  get_complement <- function(base) {
    dplyr::case_when(
      base == "A" ~ "T",
      base == "T" ~ "A",
      base == "C" ~ "G",
      base == "G" ~ "C",
      TRUE ~ base  # For non-standard bases, keep as is
    )
  }

  # Extract necessary columns as vectors
  CHR <- dat[[chr_col]]
  BP <- dat[[pos_col]]
  REF <- dat[[a1_col]]; REF_comp <- get_complement(REF)
  ALT <- dat[[a2_col]]; ALT_comp <- get_complement(ALT)
  # Construct indices
  index1 <- paste0(CHR, ":", BP, ":", REF, ":", ALT)       # Original
  index2 <- paste0(CHR, ":", BP, ":", ALT, ":", REF)       # Reversed
  index3 <- paste0(CHR, ":", BP, ":", REF_comp, ":", ALT_comp) # Complement
  index4 <- paste0(CHR, ":", BP, ":", ALT_comp, ":", REF_comp)  # Complement reversed

  # Match indices & Get rsID
  match_idx1 <- match(index1, local_ref$ID)
  match_idx2 <- match(index2, local_ref$ID)
  match_idx3 <- match(index3, local_ref$ID)
  match_idx4 <- match(index4, local_ref$ID)
  dat <- dat %>%
    dplyr::mutate(
      rsID = dplyr::case_when(
        !is.na(match_idx1) ~ local_ref$SNP[match_idx1],
        !is.na(match_idx2) ~ local_ref$SNP[match_idx2],
        !is.na(match_idx3) ~ local_ref$SNP[match_idx3],
        !is.na(match_idx4) ~ local_ref$SNP[match_idx4],
        TRUE ~ NA_character_
      )
    )

  # Remove auxiliary columns and reorder
  dat <- dat %>% dplyr::select(rsID, dplyr::everything())

  # Message about unmatched rsIDs
  na_count <- sum(is.na(dat$rsID))
  if (na_count > 0) {
    na_percent <- (na_count / nrow(dat)) * 100
    leo_message(paste0(na_count, " rsID(s) could not be matched and are set to NA. (", round(na_percent, 2), "%)"))
  }

  return(dat)
}


#' Convert rsID to CHR & BP
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
    dplyr::rename(CHR = seqnames, POS = pos, SNP = RefSNP_id) %>%
    unique()

  # Merge the SNP info with the original dataset, keeping all original columns
  merged_data <- merge(dat, snp_info, by.x = snp_col, by.y = "SNP", all.x = TRUE) %>%
    dplyr::select(SNP, CHR, POS, everything())

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
#' TODO: Merge all gene annotation function into one simple command.
#' This function maps gene symbols to their genomic positions (chromosome,
#' start, end, strand) using the specified method and genome assembly.
#'
#' The function supports two methods:
#' - `"bioconductor"`: uses Bioconductor packages. See [map_gene_to_chrbp_using_TxDb()].
#' - `"gtf"`: uses a GTF file. See [map_gene_to_chrbp_using_gtf()].
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
#' leo_map_GtoCP(genes = c("TP53", "BRCA1", "EGFR"), method = "bioconductor", genome = "hg19")
#'
#' # Using Bioconductor method with a data frame
#' leo_map_GtoCP(genes = data.frame(gene_name = c("TP53", "BRCA1", "EGFR"), value = c(1.2, 3.4, 5.6)),
#'               gene_col = "gene_name", method = "bioconductor", genome = "hg19")
#'
#' # Using GTF method with a character vector of gene symbols
#' leo_map_GtoCP(genes = c("TP53", "BRCA1", "EGFR"), method = "gtf", genome = "hg38")
#'
#' # Using GTF method with a data frame
#' leo_map_GtoCP(genes = data.frame(gene_name = c("TP53", "BRCA1", "EGFR"), value = c(1.2, 3.4, 5.6)),
#'               gene_col = "gene_name", method = "gtf", genome = "hg38", download_dir = "~/project/ref/gtf")
#' }
#' @export
leo_map_GtoCP <- function(genes,
                          gene_col = NULL,
                          method = c("bioconductor", "gtf"),
                          genome = c("hg19", "hg38"),
                          ...) {
  method <- match.arg(method)
  genome <- match.arg(genome)

  if (method == "bioconductor") {
    result <- map_gene_to_chrbp_using_TxDb(genes = genes, gene_col = gene_col, genome = genome)
  } else if (method == "gtf") {
    result <- map_gene_to_chrbp_using_gtf(genes = genes, gene_col = gene_col, genome = genome, ...)
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
#'  - For hg19, it needs `"TxDb.Hsapiens.UCSC.hg19.knownGene"` Bioconductor package
#'  - For hg38, it needs `"TxDb.Hsapiens.UCSC.hg38.knownGene"` Bioconductor package
#' @return A data frame with mapped results.
#' @importFrom dplyr %>% mutate left_join select pull
#' @importFrom rlang sym
#' @importFrom AnnotationDbi mapIds
#' @importFrom GenomicFeatures genes
#' @importFrom tidyr drop_na
#' @examples
#' \dontrun{
#' gene_symbols <- c("TP53", "BRCA1", "EGFR") # Example with gene symbol vector
#' map_gene_to_chrbp_using_TxDb(genes = gene_symbols, genome = "hg19")
#'
#' gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInformation = c(1,2,3)) # Example with data frame input
#' map_gene_to_chrbp_using_TxDb(genes = gene_symbols_df, gene_col = "GeneName" , genome = "hg19")
#' }
#' @export
map_gene_to_chrbp_using_TxDb <- function(genes, gene_col = NULL, genome = c("hg19", "hg38")) {
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
    input_df <- genes %>% tidyr::drop_na(!!gene_col)
    gene_symbols <- input_df %>% dplyr::pull(!!gene_col_sym)
  } else if (is.vector(genes)) {
    input_df <- data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE) %>%
      tidyr::drop_na(gene_symbol)
    gene_symbols <- input_df %>% pull(gene_symbol)
  } else {
    stop("Input 'genes' must be a vector or a data frame.")
  }
  unique_gene_symbols <- unique(gene_symbols) # Use unique gene symbols for mapping
  message(sprintf("Total input genes: %d; Unique genes: %d", length(gene_symbols), length(unique_gene_symbols)))

  # Map gene symbols to Entrez IDs
  entrez_ids <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                          keys = unique_gene_symbols,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals = "first") )
  # in case there are NA


  # Merge & Mapping
  mapping_df <- data.frame(
    gene_symbol = unique_gene_symbols,
    entrez_id = entrez_ids,
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(entrez_id))
  mapping_df <- mapping_df %>% dplyr::filter(!is.na(entrez_id))

  # Retrieve genomic positions using Entrez IDs
  suppressMessages(
    gene_locations <- GenomicFeatures::genes(txdb, columns = c("gene_id"),
                                             single.strand.genes.only = TRUE)
  )
  gene_locations_df <- as.data.frame(gene_locations) %>%
    dplyr::mutate(gene_id = as.character(gene_id))

  result_df <- mapping_df %>%
    dplyr::left_join(gene_locations_df, by = c("entrez_id" = "gene_id")) %>%
    dplyr::select(gene_symbol, chr = seqnames, bp_start = start, bp_end = end, strand) %>%
    dplyr::mutate(chr = gsub("^chr", "", chr)) %>%
    dplyr::distinct()

  # Group by gene_symbol to get unique positions
  standard_chr <- c(as.character(1:22), "X", "Y", "MT")
  result_df <- result_df %>%
    dplyr::filter(chr %in% standard_chr, !is.na(bp_start), !is.na(bp_end))

  # Inform about genes that could not be mapped
  unmapped_genes <- setdiff(unique_gene_symbols, result_df$gene_symbol)
  if (length(unmapped_genes) > 0) {
    message(paste0("⬇ A total of ", length(unmapped_genes), " could not be mapped ⬇"))
    message("The following genes could not be mapped: ", paste(unmapped_genes, collapse = ", "))
  }

  # Now merge the mapping back to the original input data
  final_df <- input_df %>%
    dplyr::left_join(result_df, by = stats::setNames("gene_symbol", gene_col))

  return(final_df)
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
#' map_gene_to_chrbp_using_gtf(genes = gene_symbols, genome = "hg38")
#'
#' gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInformation = c(1,2,3))
#' map_gene_to_chrbp_using_gtf(genes = gene_symbols_df, gene_col = "GeneName" , genome = "hg19")
#' }
#' @export
map_gene_to_chrbp_using_gtf <- function(genes, gene_col = NULL, genome = c("hg19", "hg38"), gtf_file = NULL, download_dir = "~/project/ref/gtf") {
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
  unique_gene_symbols <- unique(gene_symbols) # Use unique gene symbols for mapping
  message(sprintf("Total input genes: %d; Unique genes: %d", length(gene_symbols), length(unique_gene_symbols)))

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
  gene_locations_df <- gene_locations_df %>% # gene_locations_df <- as.data.frame(gene_locations) # (debug)
    dplyr::left_join(gene_symbols_map, by = c("gene_id")) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      # Extract numeric chromosome numbers
      chr = dplyr::case_when(
        grepl("^NC_", seqnames) ~ sub("^NC_0{0,}0*([0-9]+)\\..*", "\\1", seqnames),
        TRUE ~ NA_character_
      ),
      # Replace numeric codes for X and Y chromosomes
      chr = dplyr::case_when(
        chr == "23" ~ "X",
        chr == "24" ~ "Y",
        chr == "12920" ~ "MT",
        TRUE ~ chr
      )
    )

  # Filter out non-standard chromosomes and NA positions
  standard_chr <- c(as.character(1:22), "X", "Y", "MT")
  gene_locations_df <- gene_locations_df %>%
    dplyr::filter(chr %in% standard_chr, !is.na(start), !is.na(end))

  # Filter for unique gene symbols and group
  result_df <- gene_locations_df %>%
    dplyr::filter(gene_name %in% unique_gene_symbols) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      chr = chr[1],
      bp_start = min(start),
      bp_end = max(end),
      strand = strand[1]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(gene_symbol = gene_name)

  # Inform about genes that could not be mapped
  mapped_genes <- result_df$gene_symbol
  unmapped_genes <- setdiff(unique_gene_symbols, mapped_genes)
  if (length(unmapped_genes) > 0) {
    message(paste0("⬇ A total of ", length(unmapped_genes), " could not be mapped ⬇"))
    message("The following genes could not be mapped: ", paste(unmapped_genes, collapse = ", "))
  }

  # Merge back with original input data
  final_df <- input_df %>%
    dplyr::left_join(result_df, by = stats::setNames("gene_symbol", gene_col))

  return(final_df)
}


#' Map Gene Symbols to Genomic Positions Using biomaRt
#'
#' This function queries gene symbols for their genomic positions (chromosome, start, end, strand)
#' using Ensembl's biomaRt for the specified genome assembly ("hg19" or "hg38").
#'
#' If the input is a data frame, the function retains all existing columns and
#' adds new columns with the mapping results.
#'
#' @param genes A character vector of gene symbols to query, or a data frame containing gene symbols.
#' @param gene_col The column name of gene symbols if `genes` is a data frame.
#' @param genome The genome assembly to use: "hg19" or "hg38".
#'
#' @return A data frame with the original data and new columns: chr, bp_start, bp_end, strand, gene_symbol.
#' @importFrom biomaRt useMart getBM
#' @importFrom dplyr left_join mutate filter
#' @examples
#' \dontrun{
#' # Query location of TP53, BRCA1, and EGFR genes
#' gene_symbols <- c("TP53", "BRCA1", "EGFR")
#' map_gene_to_chrbp_using_biomart(genes = gene_symbols, genome = "hg19")
#'
#' # Using a data frame with gene symbols
#' gene_symbols_df <- data.frame(GeneName = gene_symbols, OtherInfo = c(1, 2, 3))
#' map_gene_to_chrbp_using_biomart(genes = gene_symbols_df, gene_col = "GeneName", genome = "hg19")
#' }
#' @export
map_gene_to_chrbp_using_biomart <- function(genes, gene_col = NULL, genome = c("hg19", "hg38")) {
  # Check and set genome argument
  genome <- match.arg(genome)
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Install it using BiocManager::install('biomaRt').")
  }

  # Connect to Ensembl based on genome assembly
  ensembl <- if (genome == "hg19") {
    biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl",
                     host = "https://grch37.ensembl.org")
  } else {
    biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl",
                     host = "https://www.ensembl.org")
  }

  # Process input genes
  if (is.data.frame(genes)) {
    if (is.null(gene_col) || !(gene_col %in% colnames(genes))) {
      stop("Please specify a valid 'gene_col' that exists in the data frame.")
    }
    gene_symbols <- genes[[gene_col]]
    input_df <- genes
  } else if (is.character(genes)) {
    gene_symbols <- genes
    input_df <- data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE)
    gene_col <- "gene_symbol"
  } else {
    stop("Input 'genes' must be a character vector or a data frame.")
  }

  # Message about the number of input genes
  total_genes <- length(gene_symbols)
  unique_genes <- length(unique(gene_symbols))
  message(sprintf("Total input genes: %d, Unique genes: %d", total_genes, unique_genes))

  # Query gene information by symbol
  unique_gene_symbols <- unique(gene_symbols)
  gene_info <- biomaRt::getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "strand",
                   "hgnc_symbol"), # "transcription_start_site" is also a attributes, but it returns too many results.
    filters = "hgnc_symbol",
    values = unique_gene_symbols,
    mart = ensembl
  )
  # Rename and re-format
  colnames(gene_info) <- c("chr", "bp_start", "bp_end", "strand", "gene_symbol")
  gene_info$gene_symbol <- as.character(gene_info$gene_symbol)

  # Filter out non-standard chromosomes (only keep 1-22, X, Y, MT)
  # standard_chr <- c(as.character(1:22), "X", "Y", "MT")
  # gene_info <- gene_info %>% dplyr::filter(chr %in% standard_chr) %>% dplyr::filter(!is.na(chr))

  # Merge back with original input to ensure consistency
  result_df <- input_df %>%
    dplyr::left_join(gene_info, by = stats::setNames("gene_symbol", gene_col))

  # Convert strand information to "+" and "-"
  result_df <- result_df %>%
    dplyr::mutate(strand = ifelse(strand == 1, "+", "-"))

  return(result_df)
}

#' Map Ensembl Gene IDs to Genomic Positions using biomaRt
#'
#' This function uses biomaRt to retrieve genomic positions (chr, bp_start, bp_end,
#' and strand) from the Ensembl database based on Ensembl gene IDs.
#'
#' @param ensembl_ids A character vector of Ensembl gene IDs, or a data frame containing this information.
#' @param ensembl_col If `ensembl_ids` is a data frame, specify the column name containing the Ensembl gene IDs.
#' @param verbose Logical indicating whether to print the unmapped information.
#' @param genome The genome version to use: "hg19" or "hg38".
#'
#' @return A data frame containing genomic position information, including ensembl_gene_id, chr, bp_start, bp_end, strand.
#' @importFrom biomaRt useMart getBM
#' @importFrom dplyr %>% left_join select mutate filter
#' @examples
#' \dontrun{
#' # Query using Ensembl gene IDs
#' ensembl_ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648")
#' map_ensg_to_chrbp_using_biomaRt(ensembl_ids = ensembl_ids, genome = "hg19")
#'
#' # Use a data frame as input
#' gene_df <- data.frame(ensembl_id = ensembl_ids, value = c(1.2, 3.4, 5.6))
#' map_ensg_to_chrbp_using_biomaRt(ensembl_ids = gene_df, ensembl_col = "ensembl_id", genome = "hg19")
#' }
#' @export
map_ensg_to_chrbp_using_biomaRt <- function(ensembl_ids, ensembl_col = NULL, genome = c("hg19", "hg38"), verbose = F) {
  # Load required package
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Please install it using BiocManager::install('biomaRt').")
  }

  genome <- match.arg(genome)

  # Select the appropriate Ensembl database
  if (genome == "hg19") {
    ensembl <- biomaRt::useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl",
      host = "https://grch37.ensembl.org"
    )
  } else if (genome == "hg38") {
    ensembl <- biomaRt::useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl",
      host = "https://www.ensembl.org"
    )
  }

  # Process the input Ensembl IDs
  if (is.data.frame(ensembl_ids)) {
    if (is.null(ensembl_col) || !(ensembl_col %in% colnames(ensembl_ids))) {
      stop("Please specify a valid 'ensembl_col', which is a column in the data frame.")
    }
    ensg_ids <- ensembl_ids[[ensembl_col]]
    input_df <- ensembl_ids
  } else if (is.vector(ensembl_ids)) {
    ensg_ids <- ensembl_ids
    input_df <- data.frame(ensembl_gene_id = ensg_ids, stringsAsFactors = FALSE)
    ensembl_col <- "ensembl_gene_id"
  } else {
    stop("Input 'ensembl_ids' must be a character vector or a data frame containing Ensembl gene IDs.")
  }

  unique_ensg_ids <- unique(ensg_ids)
  message(sprintf("Total number of input Ensembl gene IDs: %d; Number of unique gene IDs: %d",
                  length(ensg_ids), length(unique_ensg_ids)))

  # Remove version numbers from Ensembl gene IDs
  ensg_ids_clean <- sub("\\..*", "", ensg_ids)
  input_df[[ensembl_col]] <- ensg_ids_clean
  unique_ensg_ids <- unique(ensg_ids_clean)

  # Define attributes to retrieve
  attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position",
                  "strand")

  # Retrieve gene information from Ensembl
  gene_info <- biomaRt::getBM(
    attributes = attributes,
    filters = "ensembl_gene_id",
    values = unique_ensg_ids,
    mart = ensembl
  )

  # Process strand information
  gene_info <- gene_info %>% dplyr::mutate(strand = ifelse(strand == 1, "+", "-"))
  colnames(gene_info) <- c("ensembl_gene_id", "chr", "bp_start", "bp_end", "strand")

  # Merge back to the input data frame
  result_df <- input_df %>% dplyr::left_join(gene_info, by = stats::setNames("ensembl_gene_id", ensembl_col))

  # Filter out non-standard chromosomes (optional)
  # standard_chr <- c(1:22, "X", "Y", "MT")
  # result_df <- result_df %>% dplyr::filter(chr %in% standard_chr)

  # Report unmapped Ensembl IDs
  if (verbose) {
    unmapped_ensg <- result_df %>%
      dplyr::filter(is.na(chr) | is.na(bp_start) | is.na(bp_end) | is.na(strand)) %>%
      dplyr::pull(!!rlang::sym(ensembl_col))
    if (length(unmapped_ensg) > 0) {
      message("The following Ensembl gene IDs could not be mapped to genomic positions:")
      print(unmapped_ensg)
    }
  }

  return(result_df)
}

#' Map Ensembl Gene IDs to Gene Symbols using biomaRt
#'
#' This function uses biomaRt to retrieve gene symbols based on Ensembl gene IDs.
#'
#' @param ensembl_ids A character vector of Ensembl gene IDs, or a data frame containing this information.
#' @param ensembl_col If `ensembl_ids` is a data frame, specify the column name containing the Ensembl gene IDs.
#' @param type How to handle cases where one Ensembl ID maps to multiple gene symbols. Options are "combine" (default) to combine them with a separator or "first" to only use the first symbol.
#' @param sep The separator to deal with one ensg mapped to multiple gene. Default with "/".
#' @param verbose Logical indicating whether to print the unmapped information.
#' @param genome The genome version to use: "hg19" or "hg38".
#'
#' @return A data frame containing Ensembl gene IDs and corresponding gene symbols.
#' @importFrom biomaRt useMart getBM
#' @importFrom dplyr %>% left_join mutate filter group_by summarise
#' @examples
#' \dontrun{
#' # Query using Ensembl gene IDs
#' ensembl_ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648")
#' map_ensg_to_gene_using_biomaRt(ensembl_ids = ensembl_ids, genome = "hg19")
#'
#' # Use a data frame as input
#' gene_df <- data.frame(ensembl_id = ensembl_ids, value = c(1.2, 3.4, 5.6))
#' map_ensg_to_gene_using_biomaRt(ensembl_ids = gene_df, ensembl_col = "ensembl_id", genome = "hg19")
#' }
#' @export
map_ensg_to_gene_using_biomaRt <- function(ensembl_ids, ensembl_col = NULL,
  type = c("combine", "first"), sep = "/", genome = c("hg19", "hg38"), verbose = F) {

  # Load required package
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Please install it using BiocManager::install('biomaRt').")
  }

  genome <- match.arg(genome)
  type <- match.arg(type)

  # Select the appropriate Ensembl database
  if (genome == "hg19") {
    ensembl <- biomaRt::useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl",
      host = "https://grch37.ensembl.org"
    )
  } else if (genome == "hg38") {
    ensembl <- biomaRt::useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl",
      host = "https://www.ensembl.org"
    )
  }

  # Process the input Ensembl IDs
  if (is.data.frame(ensembl_ids)) {
    if (is.null(ensembl_col) || !(ensembl_col %in% colnames(ensembl_ids))) {
      stop("Please specify a valid 'ensembl_col', which is a column in the data frame.")
    }
    ensg_ids <- ensembl_ids[[ensembl_col]]
    input_df <- ensembl_ids
  } else if (is.vector(ensembl_ids)) {
    ensg_ids <- ensembl_ids
    input_df <- data.frame(ensembl_gene_id = ensg_ids, stringsAsFactors = FALSE)
    ensembl_col <- "ensembl_gene_id"
  } else {
    stop("Input 'ensembl_ids' must be a character vector or a data frame containing Ensembl gene IDs.")
  }

  unique_ensg_ids <- unique(ensg_ids)
  message(sprintf("Total number of input Ensembl gene IDs: %d; Number of unique gene IDs: %d",
                  length(ensg_ids), length(unique_ensg_ids)))

  # Remove version numbers from Ensembl gene IDs
  ensg_ids_clean <- sub("\\..*", "", ensg_ids)
  input_df[[ensembl_col]] <- ensg_ids_clean
  unique_ensg_ids <- unique(ensg_ids_clean)

  # Define attributes to retrieve
  attributes <- c("ensembl_gene_id", "hgnc_symbol")

  # Retrieve gene information from Ensembl
  gene_info <- biomaRt::getBM(
    attributes = attributes,
    filters = "ensembl_gene_id",
    values = unique_ensg_ids,
    mart = ensembl
  )

  # Process based on 'type' parameter
  if (type == "combine") {
    gene_info <- gene_info %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(hgnc_symbol = paste(unique(hgnc_symbol), collapse = sep)) %>%
      dplyr::ungroup()
  } else if (type == "first") {
    gene_info <- gene_info %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(hgnc_symbol = first(hgnc_symbol)) %>%
      dplyr::ungroup()
  }

  # Merge back to the input data frame
  result_df <- input_df %>% dplyr::left_join(gene_info, by = stats::setNames("ensembl_gene_id", ensembl_col))
  result_df <- result_df %>% mutate(hgnc_symbol = ifelse(hgnc_symbol=="", NA, hgnc_symbol))

  # Report unmapped Ensembl IDs
  if (verbose) {
    unmapped_ensg <- result_df %>%
      dplyr::filter(is.na(hgnc_symbol)) %>%
      dplyr::pull(!!rlang::sym(ensembl_col))
    if (length(unmapped_ensg) > 0) {
      message("The following Ensembl gene IDs could not be mapped to gene symbols:")
      print(unmapped_ensg)
    }
  }

  return(result_df)
}



#' Map Ensembl IDs to Gene Symbols using `org.Hs.eg.db`
#'
#' This function maps Ensembl IDs to their corresponding gene symbols using the `org.Hs.eg.db` package.
#'
#' @param ensembl_ids A character vector of Ensembl IDs or a data frame containing Ensembl IDs.
#' @param ensembl_col The column name of Ensembl IDs if `ensembl_ids` is a data frame.
#' @return A data frame with two columns: the input Ensembl IDs and the corresponding gene symbols.
#' @importFrom dplyr %>% left_join pull
#' @importFrom rlang sym
#' @importFrom AnnotationDbi mapIds
#' @examples
#' \dontrun{
#' ensembl_ids <- c("ENSG00000141510.1", "ENSG00000012048", "ENSG00000146648") # Example with Ensembl ID vector
#' map_ensg_to_gene_using_org.Hs.eg.db(ensembl_ids = ensembl_ids)
#'
#' ensembl_ids_df <- data.frame(EnsemblID = ensembl_ids, OtherInformation = c(1,2,3)) # Example with data frame input
#' map_ensg_to_gene_using_org.Hs.eg.db(ensembl_ids = ensembl_ids_df, ensembl_col = "EnsemblID")
#' }
#' @export
map_ensg_to_gene_using_org.Hs.eg.db <- function(ensembl_ids, ensembl_col = NULL) {
  # Load required package
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("Package 'org.Hs.eg.db' is required. You can install it using BiocManager::install('org.Hs.eg.db').")
    stop("Package 'org.Hs.eg.db' is required.")
  }

  # Process input Ensembl IDs
  if (is.data.frame(ensembl_ids)) {
    if (is.null(ensembl_col) || !(ensembl_col %in% colnames(ensembl_ids))) {
      stop("Please specify a valid 'ensembl_col' that exists in the data frame.")
    }
    ensembl_col_sym <- rlang::sym(ensembl_col)
    ensembl_id_values <- ensembl_ids %>% dplyr::pull(!!ensembl_col_sym)
    input_df <- ensembl_ids
  } else if (is.vector(ensembl_ids)) {
    ensembl_id_values <- ensembl_ids
    input_df <- data.frame(ensembl_id = ensembl_id_values, stringsAsFactors = FALSE)
    ensembl_col <- "ensembl_id"
  } else {
    stop("Input 'ensembl_ids' must be a vector or a data frame.")
  }

  # Remove version numbers from Ensembl IDs
  ensembl_id_clean <- sub("\\..*", "", ensembl_id_values)
  input_df[[ensembl_col]] <- ensembl_id_clean

  unique_ensembl_ids <- unique(ensembl_id_clean)
  message(sprintf("Total input Ensembl IDs: %d; Unique Ensembl IDs: %d",
                  length(ensembl_id_clean), length(unique_ensembl_ids)))

  # Map Ensembl IDs to Gene Symbols
  gene_symbols <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique_ensembl_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )

  # Create mapping data frame
  mapping_df <- data.frame(
    ensembl_id = unique_ensembl_ids,
    gene_symbol = gene_symbols[unique_ensembl_ids],
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(gene_symbol))

  # Merge back with input data
  final_df <- dplyr::left_join(input_df, mapping_df, by = stats::setNames("ensembl_id", ensembl_col))

  return(final_df)
}

#' Map Genes to Their TSS Positions
#'
#' This function maps gene symbols to their transcription start sites (TSS) positions
#' using hg19 or hg38 genome assembly.
#'
#' @param genes A character vector of gene symbols or a data frame containing gene symbols.
#' @param gene_col The column name of gene symbols if `genes` is a data frame.
#' @param genome The genome assembly to use: `"hg19"` or `"hg38"`.
#' @param gtf_file The path to a GTF file. If `NULL`, the function will download the appropriate GTF file.
#' @param download_dir The path where you want to store the downloaded gtf file.
#' @param ... Pass to map_gene_to_chrbp_using_gtf. See [map_gene_to_chrbp_using_gtf()]
#'
#' @return A data frame with gene symbols and their TSS positions
#' @importFrom dplyr %>% select mutate case_when
#' @examples
#' \dontrun{
#' genes <- c("TP53", "BRCA1", "EGFR")
#' map_gene_to_tss_using_gtf(genes = genes, genome = "hg38")
#'
#' genes_df <- data.frame(GeneName = genes, OtherInfo = c(1,2,3))
#' map_gene_to_tss_using_gtf(genes = genes_df, gene_col = "GeneName", genome = "hg19")
#' }
#' @export
map_gene_to_tss_using_gtf <- function(genes, gene_col = NULL, genome = c("hg19", "hg38"),
                                      gtf_file = NULL, download_dir = "~/project/ref/gtf", ...) {

  # Get full mapping results using existing function
  mapping_results <- map_gene_to_chrbp_using_gtf(
    genes = genes,
    gene_col = gene_col,
    genome = genome,
    gtf_file = gtf_file,
    download_dir = download_dir,
    ...
  )

  # Calculate TSS based on strand, handling NA cases
  final_df <- mapping_results %>%
    dplyr::mutate(
      tss = dplyr::case_when(
        is.na(strand) ~ NA_integer_,
        strand == "+" ~ bp_start,  # TSS is at start for forward strand
        strand == "-" ~ bp_end,    # TSS is at end for reverse strand
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::select(-bp_start, -bp_end)  # Remove the original start/end positions

  # Report genes without strand information
  genes_no_strand <- final_df %>%
    dplyr::filter(!is.na(chr), is.na(strand)) %>%
    dplyr::pull(!!rlang::sym(ifelse(is.null(gene_col), "gene_symbol", gene_col)))

  if (length(genes_no_strand) > 0) {
    message("\n⚠️ The following genes were mapped but lack strand information (TSS could not be determined):")
    message(paste(genes_no_strand, collapse = ", "))
  }

  return(final_df)
}

#' Map Ensembl Gene IDs to TSS Using biomaRt
#'
#' This function maps Ensembl gene IDs to their transcription start sites (TSS).
#' (Note that this is inferred using the gene start and end information based on strand)
#' It uses biomaRt for the specified genome assembly ("hg19" or "hg38").
#'
#' @param ensembl_ids A character vector of Ensembl gene IDs, or a data frame containing this information.
#' @param ensembl_col If `ensembl_ids` is a data frame, specify the column name containing the Ensembl gene IDs.
#' @param genome The genome version to use: "hg19" or "hg38".
#' @param ... pass to `map_ensg_to_chrbp_using_biomaRt`. See [map_ensg_to_chrbp_using_biomaRt()]
#'
#' @return A data frame with Ensembl gene IDs and their TSS positions.
#' @importFrom dplyr %>% mutate case_when select
#' @examples
#' \dontrun{
#' ensembl_ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648")
#' map_ensg_to_tss_using_biomaRt(ensembl_ids = ensembl_ids, genome = "hg19")
#'
#' ensembl_ids_df <- data.frame(EnsemblID = ensembl_ids, OtherInfo = c(1, 2, 3))
#' map_ensg_to_tss_using_biomaRt(ensembl_ids = ensembl_ids_df, ensembl_col = "EnsemblID", genome = "hg38")
#' }
#' @export
map_ensg_to_tss_using_biomaRt <- function(ensembl_ids, ensembl_col = NULL, genome = c("hg19", "hg38"), ...) {
  mapping_results <- map_ensg_to_chrbp_using_biomaRt(
    ensembl_ids = ensembl_ids,
    ensembl_col = ensembl_col,
    genome = genome,
    ...
  )

  # Calculate TSS based on strand
  final_df <- mapping_results %>%
    dplyr::mutate(
      tss = dplyr::case_when(
        strand == "+" ~ bp_start,
        strand == "-" ~ bp_end,
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::select(-bp_start, -bp_end)  # Remove original start/end positions

  return(final_df)
}

# CpG annotation
#' Annotate CpG Sites with Gene Information
#'
#' This function annotates a vector of CpG site probe IDs by retrieving corresponding gene names
#' from the Illumina 450k annotation package. It returns a data frame with the original CpG sites
#' and their associated gene names.
#'
#' @param cpg_vector A character vector of CpG site probe IDs to be annotated.
#' @param annotation_package A character string specifying the Illumina annotation package to use.
#'   It can be one of the following:
#'   - "IlluminaHumanMethylation450kanno.ilmn12.hg19" (default)
#'   - "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
#' @return A data frame with two columns:
#'   - `CpG_Site`: The original CpG site probe IDs.
#'   - `Gene`: The associated gene names. `NA` if no gene annotation is found.
#' @importFrom minfi getAnnotation
#' @importFrom dplyr filter pull mutate select
#' @importFrom cli cli_alert_info cli_alert_success
#' @examples
#' \dontrun{
#' library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#' # Example CpG probe IDs
#' cpg_ids <- c("cg00000029", "cg00000108", "cg00000109")
#' # Annotate CpG sites
#' annotate_cpg_sites(cpg_ids, "IlluminaHumanMethylation450kanno.ilmn12.hg19")
#' }
#' @export
annotate_cpg_sites <- function(cpg_vector,
                               annotation_package = "IlluminaHumanMethylation450kanno.ilmn12.hg19") {
  cli::cli_alert_info("Loading annotation package '{annotation_package}'...")
  cli::cat_bullet("Please do pre-library the required package", bullet_col = "blue")
  # Check if annotation package is installed
  if (annotation_package == "IlluminaHumanMethylation450kanno.ilmn12.hg19") {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else if (annotation_package == "IlluminaHumanMethylationEPICanno.ilm10b4.hg19") {
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else {
    stop("Invalid annotation package. Please use either 'IlluminaHumanMethylation450kanno.ilmn12.hg19' or 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'.")
  }
  annotation_data_full <- minfi::getAnnotation(get(annotation_package))
  # mapping
  cli::cli_alert_info("Filtering annotation data for provided CpG sites...")
  annotation_subset <- annotation_data_full[cpg_vector, , drop = FALSE]
  gene_names <- annotation_subset$UCSC_RefGene_Name
  cli::cli_alert_success("Annotation completed successfully.")
  result_df <- data.frame(
    CpG_Site = cpg_vector,
    Gene = gene_names,
    stringsAsFactors = FALSE
  )
  return(result_df)
}
