# roe
## Introduction

[`Alevin-fry`](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate, and memory frugal quantification tool for preprocessing single-cell RNA-sequencing data. Detailed information can be found in the alevin-fry [pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2), and [paper](https://www.nature.com/articles/s41592-022-01408-3).

The `roe` package provides useful functions for preparing input files required by `alevin-fry`, which consists of

1. preparing the *splici* reference for the `USA` mode of alevin-fry, which will export a unspliced, a spliced, and an ambiguous molecule count for each gene within each cell.
2. downloading the preprocessed quantification result of the datasets that are publicly available on the [10x Genomics website](https://www.10xgenomics.com/resources/datasets).
## Installation
The `roe` package can be accessed from its [github repository](https://github.com/COMBINE-lab/roe). To install the `roe` package, start R and enter

```{r install_roe, eval=FALSE}

# make sure the dependencies are installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install dependencies
BiocManager::install(c("eisaR","BSgenome","fishpond"))

# install roe from github
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install roe
devtools::install_github("COMBINE-lab/roe")
```


## Preparing a splici index for quantification with alevin-fry

The USA mode in alevin-fry requires a special index reference, which is called the *splici* reference. The *splici* reference contains the spliced transcripts plus the intronic sequences of each gene. The `make_splici_txome()` function is designed to make the *splici* reference by taking a genome FASTA file and a gene annotation GTF file as the input. Details about the *splici* can be found in the Section S2 and S3 of the supplementary file of the [alevin-fry paper](https://www.nature.com/articles/s41592-022-01408-3). To run the function, you also need to specify the read length argument `read_length` of the experiment you are working on and the flank trimming length `flank_trim_length`. A final flanking length will be computed as the difference between the read_length and flank trimming length and will be attached to the ends of each intron to absorb the intron-exon junctional reads. This function is based on the `getFeatureRanges()` function in [`eisaR`](https://github.com/fmicompbio/eisaR) package.

Following is an example of calling the `make_splici_txome()` to make the *splici* index reference. The final flanking length is calculated as the difference between the read length and the flank_trim_length, i.e., $5-2=3$. This function allows you to add extra spliced and unspliced sequences to the *splici* index, which will be useful when some unannotated sequences, such as mitochondrial genes, are important for your experiment. 

```{r}
library(roe)
genome_path <- system.file("extdata/small_example_genome.fa", package = "roe")
gtf_path <- system.file("extdata/small_example.gtf", package = "roe")
extra_spliced = system.file("extdata/extra_spliced.txt", package = "roe")
extra_unspliced = system.file("extdata/extra_unspliced.txt", package = "roe")
filename_prefix = "transcriptome_splici"
output_dir = tempdir()
read_length = 5
flank_trim_length = 2

make_splici_txome(
  genome_path = genome_path,
  gtf_path = gtf_path,
  read_length = read_length,
  output_dir = output_dir,
  flank_trim_length = flank_trim_length,
  filename_prefix = filename_prefix,
  extra_spliced = extra_spliced,
  extra_unspliced = extra_unspliced,
  dedup_seqs = FALSE,
  no_flanking_merge = FALSE
)
grep("transcriptome_splici", dir(output_dir), value = TRUE)
```

The `make_splici_txome()` function has no returned value, but writes two files to your specified output directory `output_dir`. They are 
- A FASTA file that stores the extracted splici sequences.
- A three columns' transcript-name-to-gene-name file that stores the name of each transcript in the splici index reference, their corresponding gene name, and the splicing status (`S` for spliced and `U` for unspliced) of those transcripts.

### the *splici* index

The *splici* index of a given species consists of the transcriptome of the species, i.e., the spliced transcripts, and the intronic sequences of the species. Within a gene, if the flanked intronic sequences overlap with each other, the overlapped intronic sequences will be collapsed as a single intronic sequence to make sure each base will appear only once in the intronic sequences. For more detailed information, please check the Section S2 and S3 in the supplementary file of the [alevin-fry paper](https://www.nature.com/articles/s41592-022-01408-3).

## preprocessed 10x datasets

10x Genomics provides many publicly available single-cell RNA-sequencing datasets. These datasets are widely used in many scientific research projects. To avoid reinventing wheels, we processed these datasets using alevin-fry using a [nextflow-based alevin-fry workflow](https://github.com/COMBINE-lab/10x-requant), and provide here R functions to download and load the quantification result of these datasets into R as SingleCellExperiment ojects.

If one needs to check the available datasets or download the available datasets to a local folder,  `preprocessed_10x_data()` is here to help.

```R
# to return the dataframe of the information of available datasets
fetch_processed_quant()

# to download some available datasets to a local directory
# according to the row index in the available_dataset_df returned 
# from the previous command
fetch_processed_quant(dataset_ids = c(1,5,7),
                      output_dir = "processed_quant",
                      force = FALSE,
                      delete_tar = TRUE,
                      quiet = FALSE
                    )
```

We also provide a thin wrapper of the `fetch_processed_quant()` and [`fishpond::loadFry`](https://github.com/mikelove/fishpond/blob/master/R/alevin-loadFry.R) so that the fetched datasets can be directly load into R as SingleCellExperiment objects. When using this function, the `outputFormat` parameter can be specified for each fetched dataset separately by providing a list of valid `outputFormat`s, named by the corresponding dataset ids. Similarly, the `nonzero` parameter can also be specified for each dataset by providing a list of valid `nonzero`s, named by the corresponding dataset ids. The following example shows how to define `output_format` in different ways.  

```R
load_processed_quant(dataset_ids = c(1, 2),
        fetch_dir = "processed_quant",
        force = FALSE,
        delete_tar = TRUE,
        output_format = "scRNA",
#         output_format = list("scRNA", "scRNA"),
#         output_format = list("1" = list(counts = c("S", "A")),
#                               "2" = list(counts = c("S", "A"))
#                              ),
#         output_format = list("counts" = c("S", "A")),
        nonzero = FALSE,
        quiet = FALSE
)
}

```