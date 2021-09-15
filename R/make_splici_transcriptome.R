#' Construct splici transcriptome
#'
#' Construct the splici transcriptome for alevin-fry.
#'
#' @param gtf_path Character scalar providing the path to a gtf file.
#' @param genome_path Character scalar providing the path to a genome fasta
#'   file.
#' @param read_length Numeric scalar indicating the read length.
#' @param flank_trim_length Numeric scalar giving the flank trimming length. The
#'   final flank length is obtained by subtracting the \code{flank_trim_length}
#'   from the \code{read_length}.
#' @param output_dir Character scalar indicating the output directory, where
#'   the splici transcriptome files will be written.
#' @param file_name_prefix Character scalar giving the file name prefix of the
#'   generated output files. The derived flank length will be automatically
#'   appended to the provided prefix.
#' @param extra_spliced Character scalar providing the path to a fasta file
#'   with additional sequences to include among the spliced ones.
#' @param extra_unspliced Character scalar providing the path to a fasta file
#'   with additional sequences to include among the unspliced ones.
#' @param dedup_seqs Logical scalar indicating whether or not to remove
#'   duplicate sequences.
#' @param write_actual_flank Logical scalar indicating whether or not to write
#'   out the actual flank length (which may be shorter than the indicated one
#'   if the latter would mean going outside the gene boundaries). 
#'   This argument is currently under development. 
#'
#' @author Dongze He
#'
#' @export
#'
#' @return Nothing is returned, but the necessary files for the splici
#'   transcriptome are written to the designated \code{output_dir}.
#'
#' @examples
#' ## Put a small example here
#'
make_splici_txome <- function(gtf_path,
                              genome_path,
                              read_length,
                              flank_trim_length = 5,
                              output_dir,
                              file_name_prefix = "transcriptome_splici",
                              extra_spliced = NULL,
                              extra_unspliced = NULL,
                              dedup_seqs = FALSE
                              # ,write_actual_flank=FALSE
                              ) {
  
  suppressWarnings(.make_splici_txome(gtf_path=gtf_path,
                                      genome_path=genome_path,
                                      read_length=read_length,
                                      flank_trim_length = flank_trim_length,
                                      output_dir=output_dir,
                                      file_name_prefix = file_name_prefix,
                                      extra_spliced=extra_spliced,
                                      extra_unspliced=extra_unspliced
                                      ,dedup_seqs=dedup_seqs
                                      # ,write_actual_flank=FALSE
                                      )
  )
}

#' @importFrom eisaR getFeatureRanges getTx2Gene
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @import BSgenome
#' @importFrom stringr str_detect word
#' @importFrom GenomicFeatures makeTxDbFromGFF genes extractTranscriptSeqs
#' @importFrom BiocGenerics unlist relist start end
#' @importFrom GenomicRanges reduce trim
#' @importFrom S4Vectors mcols metadata split
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom Biobase rowMin rowMax
#' @importFrom rlang .data
#' @importFrom dplyr %>% mutate
#' @importFrom utils write.table

.make_splici_txome <- function(gtf_path,
                               genome_path,
                               read_length,
                               flank_trim_length = 5,
                               output_dir = ".",
                               file_name_prefix = "transcriptome_splici",
                               extra_spliced = NULL,
                               extra_unspliced = NULL,
                               dedup_seqs = FALSE
                               # ,write_actual_flank=FALSE
                               ) {
  ## TODO: Add this sentence somewhere in the documentation:
  # flank trim length is used to avoid marginal case when dealing with junctional reads
  
  ############################################################################
  # Preprocessing
  ############################################################################
  
  if (!dir.exists(output_dir)) {
    dir.create(file.path(output_dir), recursive = TRUE,
               showWarnings = FALSE)
  }
  # make sure flank_length makes sense
  flank_length <- read_length - flank_trim_length
  if (flank_length < 0) {
    stop("flank trim length is larger than read length!")
  }
  # make sure gtf file exists
  if (!file.exists(gtf_path)) {
    stop("The following file does not exist: \n", gtf_path)
  }
  # make sure fasta file exists
  if (!file.exists(genome_path)) {
    stop("The following file does not exist: \n", genome_path)
  }
  
  # output file names
  file_name_prefix <- paste0(file_name_prefix, "_fl", flank_length)
  out_fa <- file.path(output_dir, paste0(file_name_prefix, ".fa"))
  # out_t2g <- file.path(output_dir, paste0(file_name_prefix, "_t2g.tsv"))
  out_t2g3col <- file.path(output_dir, paste0(file_name_prefix,
                                              "_t2g_3col.tsv"))
  
  # load the genome sequence
  x <- Biostrings::readDNAStringSet(file.path(genome_path))
  # get the first word as the name
  names(x) <- stringr::word(names(x), 1)
  
  ############################################################################
  # Process gtf to get spliced and introns
  ############################################################################
  message("============     processing gtf to get spliced and introns     ============")
  grl <- suppressWarnings(eisaR::getFeatureRanges(
    gtf = file.path(gtf_path),
    featureType = c("spliced", "intron"),
    intronType = "separate",
    flankLength = 0,
    joinOverlappingIntrons = TRUE,
    verbose = TRUE
  ))
  
  ############################################################################
  # Get spliced related stuff
  ############################################################################
  
  spliced_idx <- names(grl) %in% S4Vectors::metadata(grl)$featurelist$spliced
  spliced_grl <- grl[spliced_idx]
  
  ############################################################################
  # Get reduced introns
  ############################################################################
  
  # identify all introns and convert to GRanges
  intron_idx <- names(grl) %in% S4Vectors::metadata(grl)$featurelist$intron
  intron_gr <- BiocGenerics::unlist(grl[intron_idx])
  
  intron_gr = .add_metadata(intron_gr, x = x)
  
  # add flanking length to each side
  intron_gr_flanked = intron_gr + flank_length
  
  # This functionality hasn't been complished yet, will not expose
  ## also make sure flanked introns are within chromosome boundary
  # GenomeInfoDb::seqlevels(intron_gr_flanked) <- GenomeInfoDb::seqlevels(x)
  # GenomeInfoDb::seqlengths(intron_gr_flanked) <- suppressWarnings(
  #     GenomeInfoDb::seqlengths(x)
  # )
  # intron_gr_flanked <- GenomicRanges::trim(intron_gr_flanked)
  
  # if (write_actual_flank) {
  #     # next ensure flanked introns are within gene boundary
  #     ## first load gene ranges from GTF
  #     txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")
  #     gs <- GenomicFeatures::genes(txdb)
  # 
  #     ## grab relative gene range info into a matrix
  #     grs <- cbind(gs$gene_id, BiocGenerics::start(gs), BiocGenerics::end(gs))
  #     colnames(grs) <- c("gene_id", "gene_start", "gene_end")
  #     rownames(grs) <- grs[, "gene_id"]
  # 
  #     ## grab relative introns and flanked introns range info into another matrix
  #     trs <- cbind(intron_gr_flanked$exon_id,
  #                  BiocGenerics::start(intron_gr),
  #                  BiocGenerics::end(intron_gr),
  #                  BiocGenerics::start(intron_gr_flanked),
  #                  BiocGenerics::end(intron_gr_flanked),
  #                  sapply(strsplit(intron_gr_flanked$exon_id, "-",
  #                                  fixed = TRUE), .subset, 1), "NA")
  #     colnames(trs) <- c("intron_id", "intron_start", "intron_end",
  #                        "flanked_start", "flanked_end", "gene_id")
  #     rownames(trs) <- trs[, "intron_id"]
  # 
  #     ## merge gene range info and intron range info matrix
  #     trs <- cbind(trs, grs[trs[, "gene_id"], , drop = FALSE])
  #     trs <- trs[, !(colnames(trs) %in% c("intron_id", "gene_id"))]
  #     class(trs) <- "numeric"
  # 
  #     ## correct out-of-boundary flanking initial positions and terminal positions
  #     trs[, "flanked_start"] <- Biobase::rowMax(trs[, c("flanked_start",
  #                                                       "gene_start")])
  #     trs[, "flanked_end"] <- Biobase::rowMin(trs[, c("flanked_end",
  #                                                     "gene_end")])
  # 
  #     ## calculate the actual flanking length attached to each side of each
  #     intron_flank_length <- cbind(
  #         rownames(trs),
  #         trs[, "intron_start"] - trs[, "flanked_start"] + 1,
  #         trs[, "flanked_end"] - trs[, "intron_end"] + 1
  #     )
  #     colnames(intron_flank_length) <- c("intron_name",
  #                                        "initial_flank_length",
  #                                        "terminal_flank_length")
  # 
  # 
  #     utils::write.table(
  #         intron_flank_length,
  #         file = file.path(output_dir, "intron_flank_length.tsv"),
  #         row.names = FALSE,
  #         col.names = TRUE,
  #         quote = FALSE,
  #         sep = "\t"
  #     )
  # 
  # }
  intron_gr_flanked = .add_metadata(intron_gr_flanked, x=x)
  
  # remake intron GRangesList
  intron_grl <- BiocGenerics::relist(intron_gr_flanked, lapply(
    structure(seq_along(intron_gr_flanked),
              names = intron_gr_flanked$transcript_id), function(i) i))
  
  
  ############################################################################
  # extract sequences from genome
  ############################################################################
  
  message("============extracting spliced and intron sequences from genome============")
  
  grl <- c(spliced_grl, intron_grl)
  
  # make sure introns don't out of boundary
  GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevels(x)
  GenomeInfoDb::seqlengths(grl) <- suppressWarnings(
    GenomeInfoDb::seqlengths(x)
  )
  grl <- GenomicRanges::trim(grl)
  
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = x,
    transcripts = grl
  )
  
  # If having duplicated sequences, only keep one
  if (dedup_seqs) {
    seqs <- unique(seqs)
    grl <- grl[names(seqs)]
  }
  
  # save some space
  rm(x)
  
  ############################################################################
  # process final outputs
  ############################################################################
  message("Writing outputs...")
  
  df <- eisaR::getTx2Gene(grl)
  # utils::write.table(df, out_t2g, sep = "\t", row.names = FALSE,
  #                    quote = FALSE, col.names = FALSE)
  ## TODO: Make this more robust to possible transcript names containing '-'
  df <- df %>%
    dplyr::mutate(gene_id = stringr::word(.data$gene_id, 1, sep = '-'),
                  status = ifelse(stringr::str_detect(.data$transcript_id, '-'),
                                  'U', 'S'))
  
  Biostrings::writeXStringSet(seqs, out_fa, format = "fasta")
  utils::write.table(df, out_t2g3col
                     , sep = "\t",
                     row.names = FALSE, quote = FALSE, col.names = FALSE)


  # optional: adding extra spliced and unspliced sequences from an fasta file
  if (!is.null(extra_spliced)) {
    if (!file.exists(extra_spliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa <- file(extra_spliced, open = "r")
      lns <- readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g files and fasta file
          txp_name <- gsub(">", "", ln)
          utils::write.table(matrix(c(txp_name, txp_name), nrow = 1),
                             file = out_t2g, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
          utils::write.table(matrix(c(txp_name, txp_name, "S"),
                                    nrow = 1),
                             file = out_t2g3col, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
          utils::write.table(ln, file = out_fa, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          utils::write.table(ln, file = out_fa, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
        }
      }
    }
  }
  
  if (!is.null(extra_unspliced)) {
    if (!file.exists(extra_unspliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa <- file(extra_unspliced, open="r")
      lns <- readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g file and fasta file
          txp_name = gsub(">", "", ln)
          # utils::write.table(matrix(c(txp_name,
          #                             paste0(txp_name, "-U")),
          #                           nrow = 1), file = out_t2g,
          #                    sep = "\t",
          #                    row.names = FALSE, quote = FALSE,
          #                    col.names = FALSE, append = TRUE)
          utils::write.table(matrix(c(txp_name, txp_name, "U"),
                                    nrow = 1),
                             file = out_t2g3col, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
          utils::write.table(ln, file = out_fa, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          utils::write.table(ln, file = out_fa, sep = "\t",
                             row.names = FALSE, quote = FALSE,
                             col.names = FALSE, append = TRUE)
        }
      }
    }
  }
  
  message("Done.")
}

# This function takes a GRanage object and its genome, then returns
.add_metadata <- function(intron_gr, x) {
  # group introns by gene, then collapse overlapping ranges
  intron_grl <- GenomicRanges::reduce(S4Vectors::split(intron_gr,
                                                       intron_gr$gene_id))
  
  # clean txp names and gene names
  intron_gr <- BiocGenerics::unlist(intron_grl)
  intron_gr$exon_rank <- 1L
  intron_gr$type <- "exon"
  ## TODO: Revisit this
  intron_gr$transcript_id <- stringr::word(names(intron_gr), 1, sep = '-')
  intron_gr$gene_id <- intron_gr$transcript_id
  intron_gr$transcript_id <- make.unique(paste0(intron_gr$transcript_id, "-I"), sep = '')
  intron_gr$gene_id <- paste0(intron_gr$gene_id, "-I")
  intron_gr$exon_id <- intron_gr$transcript_id
  ## --
  names(intron_gr) <- NULL
  S4Vectors::mcols(intron_gr) <-
    S4Vectors::mcols(intron_gr)[, c("exon_id", "exon_rank",
                                    "transcript_id", "gene_id", "type")]
  
  # make sure intron ranges are within chromosome boundary
  GenomeInfoDb::seqlevels(intron_gr) <- GenomeInfoDb::seqlevels(x)
  GenomeInfoDb::seqlengths(intron_gr) <- suppressWarnings(
    GenomeInfoDb::seqlengths(x)
  )
  intron_gr <- GenomicRanges::trim(intron_gr)
  
  
}
