test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

# Test splici
## We need to do the following steps to construct splici:
## 1. read in a fasta file
## 2. read in a GTF file
## 3. for each gene, get the intronic ranges from exon ranges plus flank-length
## 4. if intronic ranges are overlapping, collapse them
## 5. extract the intronic sequences from chr.
## 6. extract the spliced txp sequences from chr.

## To test this, we can define three txps
### gene1
#### txp1: three exons with one intron I-1, I-2. These two introns doesn't overlap with each other even after attaching flank length.
####       I-1 starts at very beginning, so adding flanking length will make the range invalid.
#### txp2: three exons with two introns I-3 I-4. They overlap with each other after attaching flank length
#### txp3: two exons with one intron I-5. I-5 overlaps with I-2.

### gene2 has are the same as gene1, but in minus strand.

## Using the txps, we can test the following:
### 1. Flanking length will not make the IRanges invalid
### 2. If two introns don't overlap after adding flanking length, then they become separate introns.
### 3. If two introns do overlap after adding flanking length, then they become a collapsed intron.
### 4. if two introns overlaps, then they become a collapsed intron.
### 5. On reversed complements, the above tests hold.

# Some of the tests are taken from eisaR https://bioconductor.org/packages/release/bioc/html/eisaR.html

# context("Reference generation")

library(roe)
library(BSgenome)
library(GenomicFeatures)
library(Biostrings)


gtf_path <- system.file("extdata/small_example.gtf", package = "roe")
genome_path <- system.file("extdata/small_example_genome.fa", package = "roe")
extra_spliced <- system.file("extdata/extra_spliced.txt", package = "roe")
extra_unspliced <- system.file("extdata/extra_unspliced.txt", package = "roe")
output_dir <- tempdir()
read_length <- 5
flank_trim_length <- 2

test_that("expected files are written", {

    # create simulated genome
    genome <- Biostrings::readDNAStringSet(file.path(genome_path))
    genome_revcompl <- Biostrings::reverseComplement(genome)
    chr1 <- as.character(genome[["chr1"]])
    chr2_rev <- as.character(genome_revcompl[["chr2"]])

    # run the function
    make_splici_txome(gtf_path = gtf_path,
                      genome_path = genome_path,
                      read_length = read_length,
                      flank_trim_length = flank_trim_length,
                      output_dir = output_dir,
                      extra_spliced = extra_spliced,
                      extra_unspliced = extra_unspliced,
                      dedup_seqs = FALSE,
                      quiet = TRUE
                      # ,write_actual_flank=FALSE
                      )

    # harvest the output
    splici_seqs <- readBStringSet(file.path(output_dir, "splici_fl3.fa"))
    t2g_3col <- read.csv(file.path(output_dir, "splici_fl3_t2g_3col.tsv"),header = FALSE, sep = "\t")
    # t2g = read.csv(file.path(output_dir, "splici_fl3_t2g.tsv"),header = FALSE, sep = "\t")
    
    # Define expected sequences
    ## For each gene,
    ### we have 3 spliced txps.
    #### txp1 ([1,2], [36, 45], [71,80]), with intron regions ([3,35], [46,70])
    #### txp2 ([46,55], [91, 100]), with intron ([56,90])
    #### txp3 ([121,130], [156, 160], [191,200]) with introns ([131,155],[161,190])
    ### By adding flanking length, the intron regions become
    ####  ([0,38], [43,73])
    ####  ([53,93])
    ####  ([128,158],[168,193])
    ### after collapsing, we will left three introns for each gene
    #### I1: [1,38]
    #### I2: [43,93]
    #### I3: [128,193]

    # Check extracted sequences
    # check txp1.1
    tx1.1 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr1"), 3),
                      ranges = IRanges(c(1, 36, 71), end = c(2, 45, 80)),
                      strand = Rle(strand(c("+")), 3)
                      )
    )
    expect_equal(as.character(splici_seqs[["tx1.1"]]),
                  paste(as.character(tx1.1), collapse = ""))

    # check txp1.2
    tx1.2 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr1"), 2),
                      ranges = IRanges(c(46, 91), end = c(55, 100)),
                      strand = Rle(strand(c("+")), 2)
                      )
    )
    expect_equal(as.character(splici_seqs[["tx1.2"]]),
                  paste(as.character(tx1.2), collapse = ""))

    # check txp1.3
    tx1.3 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr1"), 3),
                      ranges = IRanges(c(121, 156, 191),
                                        end = c(130, 160, 200)),
                      strand = Rle(strand(c("+")), 3)
                      )
    )
    expect_equal(as.character(splici_seqs[["tx1.3"]]),
                  paste(as.character(tx1.3), collapse = ""))

    # check txp2.1
    tx2.1 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr2"), 3),
                      ranges = IRanges(c(1, 36, 71), end = c(2, 45, 80)),
                      strand = Rle(strand(c("-")), 3)
                      )
    )
    expect_equal(as.character(splici_seqs[["tx2.1"]]),
                  paste(rev(as.character(tx2.1)), collapse = ""))

    # check txp2.2
    tx2.2 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr2"), 2),
                      ranges = IRanges(c(46, 91), end = c(55, 100)),
                      strand = Rle(strand(c("-")), 2)
                      )
    )
    expect_equal(as.character(splici_seqs[["tx2.2"]]),
                  paste(rev(as.character(tx2.2)), collapse = ""))

    # check txp2.3
    tx2.3 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr2"), 3),
                      ranges = IRanges(c(121, 156, 191),
                                        end = c(130, 160, 200)),
                      strand = Rle(strand(c("-")), 3)
                      )
    )
    expect_equal(as.character(splici_seqs[["tx2.3"]]), 
                  paste(rev(as.character(tx2.3)), collapse = ""))

    # check g1-I
    g1_I <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr1"), 1),
                      ranges = IRanges(1, end = 38),
                      strand = Rle(strand(c("+")), 1)
                      )
    )
    expect_equal(as.character(splici_seqs[["g1-I"]]),
                  as.character(g1_I))

    # check g1-I1
    g1_I1 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr1"), 1),
                      ranges = IRanges(43, end = 93),
                      strand = Rle(strand(c("+")), 1)
                      )
    )
    expect_equal(as.character(splici_seqs[["g1-I1"]]),
                  as.character(g1_I1))

    # check g1-I2
    g1_I2 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr1"), 1),
                      ranges = IRanges(128, end = 193),
                      strand = Rle(strand(c("+")), 1)
                      )
    )
    expect_equal(as.character(splici_seqs[["g1-I2"]]),
                  as.character(g1_I2))

    # check g2-I
    g2_I <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr2"), 1),
                      ranges = IRanges(1, end = 38),
                      strand = Rle(strand(c("-")), 1)
                      )
    )
    expect_equal(as.character(splici_seqs[["g2-I"]]),
                  as.character(g2_I))

    # check g2-I1
    g2_I1 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr2"), 1),
                      ranges = IRanges(43, end = 93),
                      strand = Rle(strand(c("-")), 1)
                      )
    )
    expect_equal(as.character(splici_seqs[["g2-I1"]]),
                  as.character(g2_I1))

    # check g2-I2
    g2_I2 <- getSeq(x = genome,
                    GRanges(
                      seqnames = Rle(c("chr2"), 1),
                      ranges = IRanges(128, end = 193),
                      strand = Rle(strand(c("-")), 1)
                      )
    )
    expect_equal(as.character(splici_seqs[["g2-I2"]]), as.character(g2_I2))

    # Check t2g records
    ## t2g records are sorted, so it will be
    ### g1's txps: tx1.1, tx1.2, tx1.3
    ### g2's txps:  tx2.1, tx2.2, tx2.3
    ### g1's introns:  g1-I, g1-I1, g1-I2
    ### g2's introns:  g2-I, g2-I1, g2-I2

    ## the first column of t2g file stores txp names

    expect_equal(t2g_3col$V1, c( "tx1.1", "tx1.2", "tx1.3",
                                  "tx2.1", "tx2.2", "tx2.3",
                                  "g1-I", "g1-I1", "g1-I2",
                                  "g2-I", "g2-I1", "g2-I2",
                                  "ExtraSpliced", "ExtraUnspliced"
                                  )
                )
    
    
    
    ## the second column of t2g file stores gene names plus types (either txp, introns(-I) or unspliced(-U))

    expect_equal(t2g_3col$V2, c("g1", "g1", "g1",
                                "g2", "g2", "g2",
                                "g1", "g1", "g1",
                                "g2", "g2", "g2",
                                "ExtraSpliced", "ExtraUnspliced"
                                )
                )

    ## the third column of t2g_3col file stores the splicing status of each record, either unspliced (U) or spliced (S)
    expect_equal(t2g_3col$V3, c("S", "S", "S",
                                "S", "S", "S",
                                "U", "U", "U",
                                "U", "U", "U",
                                "S", "U"
                                )
                )

})
