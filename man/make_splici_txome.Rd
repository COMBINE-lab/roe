% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_splici_txome.R
\name{make_splici_txome}
\alias{make_splici_txome}
\title{Construct splici transcriptome}
\usage{
make_splici_txome(
  genome_path,
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
)
}
\arguments{
\item{genome_path}{Character scalar providing the path to a genome fasta
file.}

\item{gtf_path}{Character scalar providing the path to a gtf file.}

\item{read_length}{Numeric scalar indicating the read length.}

\item{output_dir}{Character scalar indicating the output directory, where
the splici transcriptome files will be written.}

\item{flank_trim_length}{Numeric scalar giving the flank trimming length. The
final flank length is obtained by subtracting the \code{flank_trim_length}
from the \code{read_length}.}

\item{filename_prefix}{Character scalar giving the file name prefix of the
generated output files. The derived flank length will be automatically
appended to the provided prefix.}

\item{extra_spliced}{Character scalar providing the path to a fasta file
with additional sequences to include among the spliced ones.}

\item{extra_unspliced}{Character scalar providing the path to a fasta file
with additional sequences to include among the unspliced ones.}

\item{dedup_seqs}{Logical scalar indicating whether or not to remove
duplicate sequences.#'}

\item{no_flanking_merge}{Logical scalar indicating whether or not to 
merge overlapping introns caused by adding flanking length.}

\item{quiet}{logical whether to display no messages}
}
\value{
Nothing is returned, but the necessary files for the splici
  reference are written to the designated \code{output_dir}.
}
\description{
Construct the splici transcriptome for alevin-fry.
}
\details{
The term \emph{splici} is shorthand for spliced + intronic reference sequence.   
This reference sequence is prepared by extracting the spliced transcripts from the reference genome according to  the desired annotation (e.g. unfiltered, or filtered for certain classes of transcripts),as well as the collapsed intervals corresponding to the introns of genes. 

The intronic sequences of the \emph{splici} reference play important roles in the various kinds of experiments discussed in the manuscript of \code{alevin-fry}.
For single-cell RNA-seq data, although one typically focuses on the  fragments and UMIs arising from the (spliced) transcriptome,and only considers the spliced and ambiguous counts when performing downstream analyses, the intronic sequences act similarly to decoy sequences proposed by Srivastava 2020 et al. 
They account for fragments that might otherwise selectively-align or pseudoalign to the transcriptome with lower quality, but that, in fact, derive from some unspliced RNA transcript molecule. 
\emph{splici} improves the accuracy and false discovery of mapping, as the intronic sequences of \code{splici} act as the decoy sequence to absorb reads map to the spliced transcriptome with low quality.
For single-nucleus RNA-seq and RNA velocity analysis, the usage is straightforward, as the intronic sequences enable \code{alevin-fry} to detect the signals from unspliced, immature RNA transcripts, which are crucial in those analyses.

To construct the \emph{splici} reference, one provides a genome fasta file and a gene annotation GTF file of the  target species, as well as the read length of the experiments  to be processed. 
Here we describe the procedure to extract the  \emph{splici} sequence of one gene; the procedure will be applied for  all genes specified in the provided GTF file.   
First, the spliced transcript sequences are extracted for all isoforms of the gene. 
Next, the intronic regions of each isoform of  the gene are extracted.  
These intronic regions consist of the  full length of the unspliced transcript, subtracting out any  exonic sequence.  
If alternative \emph{splicing} positions occur within  exons, certain sub-intervals of a gene can appear as part of both  spliced transcripts and as intronic sequence. 
Additionally, each  extracted intronic interval is extended by the provided \emph{flanking length} at its starting and ending position. 
By adding this \emph{flanking  length}, the reference can account for reads that map to the junction  of an exon and an intron.  
Otherwise these reads cannot be confidently  mapped to either the unflanked intronic region or to the spliced  isoforms of the gene. 

The intronic regions of different isoforms often overlap considerably.   
To avoid repetitive indexing of shared sub-intervals, the flanked  intronic regions across all isoforms of a gene will be combined if  they share overlapping sub-regions.  
As the result, each unique  genomic sequence will appear only once in the combined intronic  regions of the gene.  
The genomic sequences of those combined  intronic regions will then be extracted as the intronic part of this gene. 
These combined intronic sequences will each be given  a unique name, and all will be added to the \emph{splici}  reference.

The process described above will be applied to all genes  defined in the provided GTF file.  
These sequences will be  combined with any custom (user-provided) sequences, for example, one can include the sequence of  mitochondrial genes to avoid spurious mapping further. 
The resulting sequences constitute the \emph{splici} reference. 
Note that sometimes one genomic region may belong to more than one gene. 
In this case, extracted sequences can appear multiple times in \code{splici} reference with different names. 
As duplicated sequences will anyway appear only once in the \code{salmon} index except if one specifically sets the \code{-{}-keepDuplicates} flag (which we have not), we did not explicitly deduplicate repeated sequences across genes when constructing the \emph{splici} reference.  
However, this function accepts an optional argument, \code{dedup_seqs} that will perform this identical sequence deduplication during reference construction. 

Currently, the \emph{splici} index depends on both the reference genome and annotation to be indexed, as well as the length of the reads that will be mapped against this index.  
The read length is used to determine an appropriate flanking size (default of read length - 5) to add to the ends of the extracted collapsed intron sequences.  
The read length used for construction need not exactly match those being mapped, and we have evaluated that the index is reasonably robust to similar read lengths and degrades gracefully as the indexed and mapped read lengths diverge.  
Further, the short-read sequencing by synthesis most commonly employed for single-cell and single-nucleus experiments comprises a small set of characteristic read lengths.  
Nonetheless, by simply constructing the index with a single flanking length that is as large as the longest reads to be considered, and by propagating the relative splice junction annotations to the mapping step, it should be possible for the index to be constructed independently of the read length parameter.  
We are currently exploring the optimal implementation of this enhancement.
}
\examples{

library(roe)

genome_path <- system.file("extdata/small_example_genome.fa", package = "roe")
gtf_path <- system.file("extdata/small_example.gtf", package = "roe")
extra_spliced = system.file("extdata/extra_spliced.txt", package = "roe")
extra_unspliced = system.file("extdata/extra_unspliced.txt", package = "roe")
output_dir = tempdir()
read_length=5
flank_trim_length = 2
filename_prefix = "splici"

# run the function
make_splici_txome(genome_path = genome_path,
                  gtf_path = gtf_path,
                  read_length = read_length,
                  output_dir = output_dir,
                  flank_trim_length = flank_trim_length,
                  filename_prefix = filename_prefix,
                  extra_spliced = extra_spliced,
                  extra_unspliced = extra_unspliced,
                  dedup_seqs = TRUE,
                  no_flanking_merge = FALSE
) 

# grep the output filenames
grep("transcriptome_splici_fl", dir(output_dir), value = TRUE)

}
\references{
Srivastava et al. (2020).
Alignment and mapping methodology influence transcript abundance estimation
\url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8}

He et al. (2021)
Alevin-fry unlocks rapid, accurate, and memory-frugal quantification of single-cell RNA-seq data
\url{https://www.nature.com/articles/s41592-022-01408-3}
}
\author{
Dongze He
}
