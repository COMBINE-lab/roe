# Test splici
## We need to do the following steps to construct splici:
## 1. read in a fasta file
## 2. read in a GTF file
## 3. for each gene, get the intronic ranges from exon ranges plus flank-length
## 4. if intronic ranges are overlapping, collapse them
## 5. extract the intronic sequences from chr.
## 6. extract the spliced txp sequences from chr.
## 7. add in extra sequences if provided.

## So, we need to test the following aspects:
### 1. the function can return correct spliced txps.
### 2. if the introns of a gene across all txps are all independent (no overlapping sequences), 
#### - if they are overlapping after attaching flanking length, the function returns them after collapsing.
#### - if they are still independent after attaching flanking length, the function returns them without collapsing.
### 3. if the introns of a gene has overlapping sequences, the function collapse them.
### 4. if an intron of one gene has overlapping sequence with an intron of another gene, these two introns should not be collapsed.
### 5. all above cases hold for their reverse complement cases.
## To do this, we can define three genes 
### gene1's introns are all independent. but some will be overlapping if attaching frank length.
#### txp1: three exons with one intron I-1, I-2. These two introns doesn't overlap with each other even after attaching flank length.
#### txp2: two exons with one intron I-3. I-3 overlaps with I-2.
#### txp3: three exons with two introns I-4 I-5. They overlap with each other after attaching flank length

### gene2 are the same as gene1, but in reverse complement.


## Genome
genome <- Biostrings::DNAStringSet(
    c(chr1 = "TTAACATTCGCTGGGGGAGATGACGAGACTAGCCGCCGCGTGGTCCTGCCGCATTATACGTGTTCAAGCGCCTACGTGGGTTGGGCAACCCGTGCCTATGGAGGCATGGACAAATTAGGTTCAACTTCAGCTACGTACGAGACCTAGAGGTAATAAGGGTATTTTACTCGGAGCATGTTTCAGTACGAACGTTAGATATC",
      chr2 = "CTATCGAAGTGGAATCTTGAAGAGCCCATCGGTTAAGGTCTCTCCAATGTCCAGCCTATTCTATGGCACGGCAGACCCGTTGTGCATCCACAGTGATAACTTACTTGGGCTCTTAATAGAGGAGTGTTGCCATTTTATCGGCTTGCACTCCAATTAGCACCAAGTGCCGTTATTGGGGTATTGCACTCATCAATAGCGTG")
)
genome_revcompl <- Biostrings::reverseComplement(genome)

## Transcript-to-gene mapping
t2g <- data.frame(
    transcript_id = c("tx1.1", "tx1.2", "tx1.3", "tx2.1", "tx2.2", "tx2.3"),
    gene_id = c("g1", "g1", "g1", "g2", "g2", "g2"),
    stringsAsFactors = FALSE
)

## Transcripts
txome <- Biostrings::DNAStringSet(
    c(`tx1.1` = paste0(substr(genome[["chr1"]], 1, 2),
                       substr(genome[["chr1"]], 36, 45),
                       substr(genome[["chr1"]], 71, 80)),
      `tx1.2` = paste0(substr(genome[["chr1"]], 46, 55),
                       substr(genome[["chr1"]], 91, 100)),
      `tx1.3` = paste0(substr(genome[["chr1"]], 121, 130),
                       substr(genome[["chr1"]], 156, 160),
                       substr(genome[["chr1"]], 191, 200)),
      `tx2.1` = paste0(substr(genome_revcompl[["chr2"]], 1, 2),
                       substr(genome_revcompl[["chr2"]], 36, 45),
                       substr(genome_revcompl[["chr2"]], 71, 80)),
      `tx2.2` = paste0(substr(genome_revcompl[["chr2"]], 46, 55),
                       substr(genome_revcompl[["chr1"]], 91, 100),
      `tx2.3` = paste0(substr(genome_revcompl[["chr2"]], 121, 130),
                       substr(genome_revcompl[["chr2"]], 156, 160),
                       substr(genome_revcompl[["chr2"]], 191, 200)))
    )
)

## GTF
gtf <- GenomicRanges::GRanges(
    seqnames = rep(c("chr1", "chr2"), c(12, 12)),
    ranges = IRanges::IRanges(
        start = c(1, 36, 71, # exons of txp1
                  46, 91, # exons of txp2
                  121, 156, 191, # exons of txp3
                  1, 46, 121, # txps
                  1, # gene1
                  1, 36, 71, # exons of txp1
                  46, 91, # exons of txp2
                  121, 156, 191, # exons of txp3
                  1, 46, 121, # txps
                  1 # gene2
                  ),
        end = c(2, 45, 80, # exons of txp1
                55, 100, # exons of txp2
                130, 160, 200, # exons of txp3
                80, 100, 200, # txps
                200, # gene1
                2, 45, 80, # exons of txp1
                55, 100, # exons of txp2
                130, 160, 200, # exons of txp3
                80, 100, 200, # txps
                200 # gene2
                )
    ),
    strand = rep(c("+", "-"), c(12, 12)),
    type = c(rep("exon", 8), rep("transcript", 3), "gene", rep("exon", 8), rep("transcript", 3), "gene"),
    gene_id = rep(c("g1", "g2"), c(12,12)),
    transcript_id = c(rep("tx1.1", 3), rep("tx1.2", 2), rep("tx1.3", 3),"tx1.1", "tx1.2", "tx1.3", NA, 
                      rep("tx2.1", 3), rep("tx2.2", 2), rep("tx2.3", 3),"tx2.1", "tx2.2", "tx2.3", NA),
    exon_id = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", NA, NA, NA, NA,
                "E9", "E10", "E11", "E12", "E13", "E14", "E15", "E16", NA, NA, NA, NA)
)

rtracklayer::export(gtf, con = "small_example.gtf", format = "gtf")
Biostrings::writeXStringSet(genome, filepath = "small_example_genome.fa")


### Finally, let's make some extra sequences
extra_spliced = c(">ExtraSpliced", "ATATATATATATATATATATATATATATATATATATATAT")
extra_unspliced = c(">ExtraUnspliced", "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")

write.table(extra_spliced, file = "extra_spliced.txt",quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
write.table(extra_unspliced, file = "extra_unspliced.txt",quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)


