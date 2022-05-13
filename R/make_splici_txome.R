#' Construct splici transcriptome
#'
#' Construct the splici transcriptome for alevin-fry.
#'
#' @param genome_path Character scalar providing the path to a genome fasta
#'   file.
#' @param gtf_path Character scalar providing the path to a gtf file.
#' @param read_length Numeric scalar indicating the read length.
#' @param flank_trim_length Numeric scalar giving the flank trimming length. The
#'   final flank length is obtained by subtracting the \code{flank_trim_length}
#'   from the \code{read_length}.
#' @param output_dir Character scalar indicating the output directory, where
#'   the splici transcriptome files will be written.
#' @param filename_prefix Character scalar giving the file name prefix of the
#'   generated output files. The derived flank length will be automatically
#'   appended to the provided prefix.
#' @param extra_spliced Character scalar providing the path to a fasta file
#'   with additional sequences to include among the spliced ones.
#' @param extra_unspliced Character scalar providing the path to a fasta file
#'   with additional sequences to include among the unspliced ones.
#' @param dedup_seqs Logical scalar indicating whether or not to remove
#'   duplicate sequences.#' 
#' @param no_flanking_merge Logical scalar indicating whether or not to 
#'   merge overlapping introns caused by adding flanking length.
#' @param quiet logical whether to display no messages
#'
#' @author Dongze He
#'
#'
#' @details
#' 
#' The term \emph{splici} is shorthand for spliced + intronic reference sequence.   
#' This reference sequence is prepared by extracting the spliced transcripts from the reference genome according to  the desired annotation (e.g. unfiltered, or filtered for certain classes of transcripts),as well as the collapsed intervals corresponding to the introns of genes. 
#' 
#' The intronic sequences of the \emph{splici} reference play important roles in the various kinds of experiments discussed in the manuscript of \code{alevin-fry}.
#' For single-cell RNA-seq data, although one typically focuses on the  fragments and UMIs arising from the (spliced) transcriptome,and only considers the spliced and ambiguous counts when performing downstream analyses, the intronic sequences act similarly to decoy sequences proposed by Srivastava 2020 et al. 
#' They account for fragments that might otherwise selectively-align or pseudoalign to the transcriptome with lower quality, but that, in fact, derive from some unspliced RNA transcript molecule. 
#' \emph{splici} improves the accuracy and false discovery of mapping, as the intronic sequences of \code{splici} act as the decoy sequence to absorb reads map to the spliced transcriptome with low quality.
#' For single-nucleus RNA-seq and RNA velocity analysis, the usage is straightforward, as the intronic sequences enable \code{alevin-fry} to detect the signals from unspliced, immature RNA transcripts, which are crucial in those analyses.
#' 
#' To construct the \emph{splici} reference, one provides a genome fasta file and a gene annotation GTF file of the  target species, as well as the read length of the experiments  to be processed. 
#' Here we describe the procedure to extract the  \emph{splici} sequence of one gene; the procedure will be applied for  all genes specified in the provided GTF file.   
#' First, the spliced transcript sequences are extracted for all isoforms of the gene. 
#' Next, the intronic regions of each isoform of  the gene are extracted.  
#' These intronic regions consist of the  full length of the unspliced transcript, subtracting out any  exonic sequence.  
#' If alternative \emph{splicing} positions occur within  exons, certain sub-intervals of a gene can appear as part of both  spliced transcripts and as intronic sequence. 
#' Additionally, each  extracted intronic interval is extended by the provided \emph{flanking length} at its starting and ending position. 
#' By adding this \emph{flanking  length}, the reference can account for reads that map to the junction  of an exon and an intron.  
#' Otherwise these reads cannot be confidently  mapped to either the unflanked intronic region or to the spliced  isoforms of the gene. 
#' 
#' The intronic regions of different isoforms often overlap considerably.   
#' To avoid repetitive indexing of shared sub-intervals, the flanked  intronic regions across all isoforms of a gene will be combined if  they share overlapping sub-regions.  
#' As the result, each unique  genomic sequence will appear only once in the combined intronic  regions of the gene.  
#' The genomic sequences of those combined  intronic regions will then be extracted as the intronic part of this gene. 
#' These combined intronic sequences will each be given  a unique name, and all will be added to the \emph{splici}  reference.
#' 
#' The process described above will be applied to all genes  defined in the provided GTF file.  
#' These sequences will be  combined with any custom (user-provided) sequences, for example, one can include the sequence of  mitochondrial genes to avoid spurious mapping further. 
#' The resulting sequences constitute the \emph{splici} reference. 
#' Note that sometimes one genomic region may belong to more than one gene. 
#' In this case, extracted sequences can appear multiple times in \code{splici} reference with different names. 
#' As duplicated sequences will anyway appear only once in the \code{salmon} index except if one specifically sets the \code{-{}-keepDuplicates} flag (which we have not), we did not explicitly deduplicate repeated sequences across genes when constructing the \emph{splici} reference.  
#' However, this function accepts an optional argument, \code{dedup_seqs} that will perform this identical sequence deduplication during reference construction. 
#' 
#' Currently, the \emph{splici} index depends on both the reference genome and annotation to be indexed, as well as the length of the reads that will be mapped against this index.  
#' The read length is used to determine an appropriate flanking size (default of read length - 5) to add to the ends of the extracted collapsed intron sequences.  
#' The read length used for construction need not exactly match those being mapped, and we have evaluated that the index is reasonably robust to similar read lengths and degrades gracefully as the indexed and mapped read lengths diverge.  
#' Further, the short-read sequencing by synthesis most commonly employed for single-cell and single-nucleus experiments comprises a small set of characteristic read lengths.  
#' Nonetheless, by simply constructing the index with a single flanking length that is as large as the longest reads to be considered, and by propagating the relative splice junction annotations to the mapping step, it should be possible for the index to be constructed independently of the read length parameter.  
#' We are currently exploring the optimal implementation of this enhancement.
#' 
#' @references
#' Srivastava et al. (2020).
#' Alignment and mapping methodology influence transcript abundance estimation
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8}
#' 
#' He et al. (2021)
#' Alevin-fry unlocks rapid, accurate, and memory-frugal quantification of single-cell RNA-seq data
#' \url{https://www.nature.com/articles/s41592-022-01408-3}
#' 
#'
#'
#' @export
#'
#' @return Nothing is returned, but the necessary files for the splici
#'   reference are written to the designated \code{output_dir}.
#'
#' @examples
#' 
#' library(roe)
#' 
#' genome_path <- system.file("extdata/small_example_genome.fa", package = "roe")
#' gtf_path <- system.file("extdata/small_example.gtf", package = "roe")
#' extra_spliced = system.file("extdata/extra_spliced.txt", package = "roe")
#' extra_unspliced = system.file("extdata/extra_unspliced.txt", package = "roe")
#' output_dir = tempdir()
#' read_length=5
#' flank_trim_length = 2
#' filename_prefix = "transcriptome_splici"
#' 
#' # run the function
#' make_splici_txome(genome_path = genome_path,
#'                   gtf_path = gtf_path,
#'                   read_length = read_length,
#'                   output_dir = output_dir,
#'                   flank_trim_length = flank_trim_length,
#'                   filename_prefix = filename_prefix,
#'                   extra_spliced = extra_spliced,
#'                   extra_unspliced = extra_unspliced,
#'                   dedup_seqs = TRUE,
#'                   no_flanking_merge = FALSE
#' ) 
#' 
#' # grep the output filenames
#' grep("transcriptome_splici_fl", dir(output_dir), value = TRUE)
#' 


make_splici_txome <- function(genome_path,
                              gtf_path,
                              read_length,
                              output_dir,
                              flank_trim_length = 5,
                              filename_prefix = "splici",
                              extra_spliced = NULL,
                              extra_unspliced = NULL,
                              dedup_seqs = FALSE,
                              no_flanking_merge = FALSE,
                              quiet = FALSE
                              ) {

  suppressWarnings(.make_splici_txome(genome_path = genome_path,
                                      gtf_path = gtf_path,
                                      read_length = read_length,
                                      flank_trim_length = flank_trim_length,
                                      output_dir = output_dir,
                                      filename_prefix = filename_prefix,
                                      extra_spliced = extra_spliced,
                                      extra_unspliced = extra_unspliced,
                                      dedup_seqs = dedup_seqs,
                                      no_flanking_merge = no_flanking_merge,
                                      quiet = FALSE
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

.make_splici_txome <- function(genome_path,
                                gtf_path,
                                read_length,
                                output_dir,
                                flank_trim_length = 5,
                                filename_prefix = "splici",
                                extra_spliced = NULL,
                                extra_unspliced = NULL,
                                dedup_seqs = FALSE,
                                no_flanking_merge = FALSE,
                                quiet = FALSE
                                ) {
  ## TODO: Add this sentence somewhere in the documentation:
  # flank trim length is used to avoid marginal case when
  # dealing with junctional reads

  ############################################################################
  # Preprocessing
  ############################################################################

  .say(quiet, "- Locating required files...")
  if (!dir.exists(output_dir)) {
    dir.create(file.path(output_dir), recursive = TRUE,
                showWarnings = FALSE)
  }
  # make sure flank_length makes sense
  flank_length <- read_length - flank_trim_length
  if (flank_length < 0) {
    stop("The flank trim length must be smaller than the read length!")
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
  filename_prefix <- paste0(filename_prefix, "_fl", flank_length)
  out_fa <- file.path(output_dir, paste0(filename_prefix, ".fa"))
  # out_t2g <- file.path(output_dir, paste0(filename_prefix, "_t2g.tsv"))
  out_t2g3col <- file.path(output_dir, paste0(filename_prefix,
                                              "_t2g_3col.tsv"))
  out_dup <- file.path(output_dir, "duplicate_entries.tsv")
  .say(quiet, "- Processing spliced transcripts and introns...")

  .say(quiet, "  - Loading the input files")
  # load the genome sequence
  x <- Biostrings::readDNAStringSet(file.path(genome_path))
  # get the first word as the name
  names(x) <- stringr::word(names(x), 1)

  ############################################################################
  # Process gtf to get spliced transcripts and introns
  ############################################################################

  suppressMessages({
    suppressWarnings({grl <- eisaR::getFeatureRanges(
      gtf = file.path(gtf_path),
      featureType = c("spliced", "intron"),
      intronType = "separate",
      flankLength = 0,
      joinOverlappingIntrons = TRUE,
      verbose = FALSE
    )})
  })

  ############################################################################
  # Get spliced related stuff
  ############################################################################

  .say(quiet, "  - Processing spliced transcripts")
  spliced_idx <- names(grl) %in% S4Vectors::metadata(grl)$featurelist$spliced
  spliced_grl <- grl[spliced_idx]

  ############################################################################
  # Get reduced introns
  ############################################################################
  .say(quiet, "  - Processing introns")

  # identify all introns and convert to GRanges
  intron_idx <- names(grl) %in% S4Vectors::metadata(grl)$featurelist$intron
  intron_gr <- BiocGenerics::unlist(grl[intron_idx])
  # pre-flanking merge
  if (no_flanking_merge) {
    intron_gr <- .add_metadata(intron_gr, x = x)
  }

  # add flanking length to each side
  intron_gr_flanked <- intron_gr + flank_length

  if (!no_flanking_merge) {
    intron_gr_flanked <- .add_metadata(intron_gr_flanked, x =x)
  }

  # remake intron GRangesList
  intron_grl <- BiocGenerics::relist(intron_gr_flanked,
                    lapply(structure(seq_along(intron_gr_flanked),
                                    names = intron_gr_flanked$transcript_id
                          ),
                          function(i) i
                    )
                )

  ############################################################################
  # extract sequences from genome
  ############################################################################
  .say(quiet, "  - Extracting sequences from the genome")

  grl <- c(spliced_grl, intron_grl)

  # make sure introns don't out of boundary
  suppressWarnings({
    GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevels(x)
    GenomeInfoDb::seqlengths(grl) <- GenomeInfoDb::seqlengths(x)
    grl <- GenomicRanges::trim(grl)
  })
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = x,
    transcripts = grl
  )

  # If having duplicated sequences, only keep one
  if (dedup_seqs) {
    seqs <- dedup_sequences(seqs)
    grl <- grl[names(seqs)]
  }

  # save some space
  rm(x)

  ############################################################################
  # process final outputs
  ############################################################################
  .say(quiet, "- Writing outputs")


  df <- eisaR::getTx2Gene(grl)
  # utils::write.table(df, out_t2g, sep = "\t", row.names = FALSE,
  #                    quote = FALSE, col.names = FALSE)
  ## TODO: Make this more robust to possible transcript names containing '-'
  df <- df %>%
    dplyr::mutate(gene_id = stringr::word(.data$gene_id, 1, sep = '-'),
                  status = ifelse(stringr::str_detect(.data$transcript_id, '-'),
                                  'U', 'S'))
  .say(quiet, "  - Writing spliced and intron sequences")

  Biostrings::writeXStringSet(seqs, out_fa, format = "fasta")
  utils::write.table(df, out_t2g3col
                      , sep = "\t",
                      row.names = FALSE, quote = FALSE, col.names = FALSE)

  # optional: adding extra spliced and unspliced sequences from an fasta file
  if (!is.null(extra_spliced)) {
    .say(quiet, "  - Appending extra spliced sequences")

    if (!file.exists(extra_spliced)) {
      warning("  - Provided extra_sequences file does not exist, ignored.")
    } else {
      fa <- file(extra_spliced, open = "r")
      lns <- readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g files and fasta file
          txp_name <- gsub(">", "", ln)
          # utils::write.table(matrix(c(txp_name, txp_name), nrow = 1),
          #                    file = out_t2g, sep = "\t",
          #                    row.names = FALSE, quote = FALSE,
          #                    col.names = FALSE, append = TRUE)
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
    .say(quiet, "  - Appending extra unspliced sequences")

    if (!file.exists(extra_unspliced)) {
      warning("  - Provided extra_sequences file does not exist, ignored.")
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
  .say(quiet, "Done")
}

# This function takes a GRanage object and its genome, then returns
.add_metadata <- function(intron_gr, x) {
  # group introns by gene, then collapse overlapping ranges
  intron_grl <- GenomicRanges::reduce(S4Vectors::split(intron_gr,
                                                        intron_gr$gene_id))
  intron_gr <- BiocGenerics::unlist(intron_grl)

  # clean txp names and gene names
  intron_gr$exon_rank <- 1L
  intron_gr$type <- "exon"
  ## TODO: Revisit this
  intron_gr$transcript_id <- stringr::word(names(intron_gr), 1, sep = '-')
  intron_gr$gene_id <- intron_gr$transcript_id
  intron_gr$transcript_id <- make.unique(paste0(intron_gr$transcript_id,
                                                "-I"),
                                        sep = "")
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
  intron_gr
}

dedup_sequences <- function(seqs, out_dup) {
  # sort seqs based on names
  # so that unique() will keep the
  # one with the smallest lex order
  seqs <- seqs[sort(names(seqs))]

  # save some work, only build list for duplicated items
  # get all non-unique seqs
  # the representative is not included
  duplicated_seqs <- seqs[duplicated(seqs)]

  # this is what we will return
  unique_seqs <- unique(seqs)
  rm(seqs)

  # we find all the names for each duplicated seq
  record_representatives <- list()
  while (length(duplicated_seqs) > 0) {
    # find the the elements that have
    # the same seq with the first element
    dup_idx <- which(duplicated_seqs == duplicated_seqs[1])
    dup_names <- names(duplicated_seqs)[dup_idx]

    # find the representative seq in the unqie_seqs
    repr_name <- names(unique_seqs)[unique_seqs == duplicated_seqs[1]]

    # record them
    record_representatives[[repr_name]] <- dup_names

    # remove them from duplicated_seqs
    duplicated_seqs <- duplicated_seqs[-dup_idx]
  }

  # dump dropped names as salmon
  out_df <- data.frame(RetainedRef = rep(names(record_representatives),
                                          sapply(record_representatives,
                                                length)
                                      ),
                        DuplicateRef= unlist(record_representatives)
            )
  write.table(out_df, file = out_dup, quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(out_df, file = "testfile", quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = TRUE)

  # return unique seqs
  return(unique_seqs)
}