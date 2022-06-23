#' Fetch preprocessed quantification result.
#'
#' Fetch alevin-fry processed quantification 
#' result of publicly available datasets.
#'
#' @param dataset_ids integer scalar or vector providing
#' the id of the available dataset(s) to be fetched.
#' @param fetch_dir path to the directory where the fetched
#' quantification results will be stored. It will be created if not exists.
#' @param force logical whether to force re-fetching the existing datasets.
#' @param delete_tar logical whether to delete the compressed datasets after
#' decompressing. If FALSE, the tar files
#'  will be stored in a folder called
#' datasets_tar under the \code{fetch_dir}.
#' @param quiet logical whether to display no messages.
#'
#' @author Dongze He
#'
#'
#' @details
#' The raw data for many single-cell and single-nucleus RNA-seq
#' experiments is publicly available. However, certain datasets
#' are used again and again, to demonstrate data processing in
#' tutorials, as benchmark datasets for novel methods (e.g. for
#' clustering, dimensionality reduction, cell type identification
#' , etc.). In particular, 10x Genomics hosts various publicly
#' available datasets generated using their technology and
#' processed via their Cell Ranger software on
#' \href{https://www.10xgenomics.com/resources/datasets}{their website}
#' for download.
#' 
#' We have created a \href{https://www.nextflow.io}{Nextflow}-based 
#' \code{alevin-fry} workflow that one can use to easily quantify
#' single-cell RNA-sequencing data in a single workflow.  The 
#' pipeline can be found \href{https://github.com/COMBINE-lab/10x-requant}{here}.
#' To test out this initial pipeline, we have begun to reprocess the
#' publicly-available datasets collected from the 10x website. We have
#' focused the initial effort on standard single-cell and single-nucleus
#' gene-expression data generated using the Chromium v2 and v3 chemistries,
#' but hope to expand the pipeline to more complex protocols soon
#' (e.g. feature barcoding experiments) and process those data as well.
#' We note that these more complex protocols can already be processed with
#' \code{alevin-fry} (see the
#' \href{https://combine-lab.github.io/alevin-fry-tutorials}{alevin-fry tutorials}),
#' but these have just not yet been incorporated into the
#' automated Nextflow-based workflow linked above.
#' 
#' Following we list the name, link and dataset id of the currently
#' available datasets whose quantification result is ready for fetch.
#' To obtain the details of these available datasets as a data frame,
#' simply run `fetch_processed_quant()` in R.
#' 
#' \enumerate{
#' \item \href{https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-controller-3-1-low-6-1-0}{500 Human PBMCs, 3' LT v3.1, Chromium Controller}
#' \item \href{https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-x-3-1-low-6-1-0}{500 Human PBMCs, 3' LT v3.1, Chromium X}
#' \item \href{https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0}{1k PBMCs from a Healthy Donor (v3 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0}{10k PBMCs from a Healthy Donor (v3 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-ht-v3-1-chromium-x-3-1-high}{10k Human PBMCs, 3' v3.1, Chromium X}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high}{10k Human PBMCs, 3' v3.1, Chromium Controller}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-single-indexed-3-1-standard-4-0-0}{10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Single Indexed}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-dual-indexed-3-1-standard-4-0-0}{10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Dual Indexed}
#' \item \href{https://www.10xgenomics.com/resources/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0}{20k Human PBMCs, 3' HT v3.1, Chromium X}
#' \item \href{https://www.10xgenomics.com/resources/datasets/pbmcs-3p_edta_sepmate-3-1-standard}{PBMCs from EDTA-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/pbmcs-3p_heparin_sepmate-3-1-standard}{PBMCs from Heparin-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/pbmcs-3p_acda_sepmate-3-1-standard}{PBMCs from ACD-A Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/pbmcs-3p_citrate_sepmate-3-1-standard}{PBMCs from Citrate-Treated Blood Collection Tubes Isolated via SepMate-Ficoll Gradient (3' v3.1 Chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/pbmcs-3p_citrate_cpt-3-1-standard}{PBMCs from Citrate-Treated Cell Preparation Tubes (3' v3.1 Chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/pbm-cs-from-a-healthy-donor-whole-transcriptome-analysis-3-1-standard-4-0-0}{PBMCs from a Healthy Donor: Whole Transcriptome Analysis}
#' \item \href{https://www.10xgenomics.com/resources/datasets/whole-blood-rbc-lysis-for-pbmcs-neutrophils-granulocytes-3-3-1-standard}{Whole Blood RBC Lysis for PBMCs and Neutrophils, Granulocytes, 3'}
#' \item \href{https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-5-3-1-standard-3-1-0}{Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 5)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-1-3-1-standard-3-1-0}{Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 1)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-chromium-connect-channel-5-3-1-standard-3-1-0}{Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Chromium Connect (channel 5)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-chromium-connect-channel-1-3-1-standard-3-1-0}{Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Chromium Connect (channel 1)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/hodgkins-lymphoma-dissociated-tumor-whole-transcriptome-analysis-3-1-standard-4-0-0}{Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis}
#' \item \href{https://www.10xgenomics.com/resources/datasets/200-sorted-cells-from-human-glioblastoma-multiforme-3-lt-v-3-1-3-1-low-6-0-0}{200 Sorted Cells from Human Glioblastoma Multiforme, 3’ LT v3.1}
#' \item \href{https://www.10xgenomics.com/resources/datasets/750-sorted-cells-from-human-invasive-ductal-carcinoma-3-lt-v-3-1-3-1-low-6-0-0}{750 Sorted Cells from Human Invasive Ductal Carcinoma, 3’ LT v3.1}
#' \item \href{https://www.10xgenomics.com/resources/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0}{2k Sorted Cells from Human Glioblastoma Multiforme, 3’ v3.1}
#' \item \href{https://www.10xgenomics.com/resources/datasets/7-5-k-sorted-cells-from-human-invasive-ductal-carcinoma-3-v-3-1-3-1-standard-6-0-0}{7.5k Sorted Cells from Human Invasive Ductal Carcinoma, 3’ v3.1}
#' \item \href{https://www.10xgenomics.com/resources/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0}{Human Glioblastoma Multiforme: 3’v3 Whole Transcriptome Analysis}
#' \item \href{https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0}{1k Brain Cells from an E18 Mouse (v3 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0}{10k Brain Cells from an E18 Mouse (v3 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0}{1k Heart Cells from an E18 mouse (v3 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0}{10k Heart Cells from an E18 mouse (v3 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-single-indexed-3-1-standard-4-0-0}{10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, Single Indexed}
#' \item \href{https://www.10xgenomics.com/resources/datasets/10-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-dual-indexed-3-1-standard-4-0-0}{10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, Dual Indexed}
#' \item \href{https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-2-chemistry-3-standard-3-0-0}{1k PBMCs from a Healthy Donor (v2 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-2-chemistry-3-standard-3-0-0}{1k Brain Cells from an E18 Mouse (v2 chemistry)}
#' \item \href{https://www.10xgenomics.com/resources/datasets/1-k-heart-cells-from-an-e-18-mouse-v-2-chemistry-3-standard-3-0-0}{1k Heart Cells from an E18 mouse (v2 chemistry)}
#' }
#' 
#' Note that because the name of datasets are too long, the stored
#' datasets are named by their id.
#'
#' @export
#'
#' @return If an empty dataset_ids is provided,
#' a data frame containing the information of
#' available datasets will be returned;
#' otherwise, a list of ProcessedQuant class objects, in which
#' each ProcessedQuant object stores the information of one
#' fetched dataset. The `quant_path` field
#' represents the path to the quantification
#' result of the fetched dataset.
#'
#' @examples
#' 
#' \dontrun{
#' library(roe)
#' # run the function
#' available_datasets = load_processed_quant()
#' fetched_quant_list = fetch_processed_quant(dataset_id = c(1, 3),
#'                                               fetch_dir = "processed_quant",
#'                                               force = FALSE,
#'                                               delete_tar = FALSE,
#'                                               quiet = FALSE)
#' 
#' print(fetched_quant_list$"1"@quant_path)
#' print(fetched_quant_list$"2"@quant_path)
#' }
#'

fetch_processed_quant <- function(dataset_ids = c(),
                                    fetch_dir = "processed_quant",
                                    force = FALSE,
                                    delete_tar = FALSE,
                                    quiet = FALSE
                                ) {

    # available_datasets = read.csv("available_datasets.tsv", sep = "\t")
    # usethis::use_data(available_datasets, internal = TRUE, force = TRUE)

    # if the user just wants the data frame, return it
    if (length(dataset_ids) == 0) {
        return(available_datasets)
    }

    .say(quiet, "Checking provided dataset ids")

    dataset_ids <- check_dataset_ids(dataset_ids)
    # check whether there is any dataset id left
    if (length(dataset_ids) == 0) {
        stop("No valid dataset id found, can not proceed")
    }

    # download the quantification tar file for each queried dataset.
    pq_list <- list()
    # folder for (temporarily) storing tar files.
    tar_dir <- file.path(fetch_dir, "quant_tar")
    dir.create(tar_dir, recursive = TRUE,
                showWarnings = FALSE)

    for (dataset_id in dataset_ids) {
        # init processed_quant
        processed_quant <- ProcessedQuant(dataset_id)

        # fetch it
        processed_quant <- fetch_quant(processed_quant,
                                    tar_dir = tar_dir,
                                    force = force,
                                    quiet = quiet)

        # decompress it
        processed_quant <- decompress_quant(processed_quant,
                                            quant_dir = fetch_dir,
                                            force = force,
                                            quiet = quiet)

        # reset tar_path if needed
        if (delete_tar) {
            processed_quant@tar_path <- character(0)
        }

        # append to list
        pq_list[[as.character(dataset_id)]] <- processed_quant
    }

    if (delete_tar) {
        .say(quiet,
            "Removing downloaded tar files in directory:\n",
            "  ", tar_dir)
        unlink(tar_dir,  recursive = TRUE, force = TRUE)
    }

    .say(quiet, "Done")
    return(pq_list)
}
