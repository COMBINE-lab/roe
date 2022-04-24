#' Query preprocessed 10x datasets
#'
#' Construct the splici transcriptome for alevin-fry.
#'
#' @param dataset_ids integer scalar or vector providing the id of the
#' dataset(s) to be downloaded and processed. If a scalar is given, the
#' returned object will be a SingleCellExperiment object; If a vector
#' is given, the returned object will be a list of SingleCellExperiment
#' objects, each for a queried dataset.
#' @param output_dir path to the output folder, will create if not exists.
#' @param force logical whether to force re-downloading the existing datasets.
#' @param delete_tar logical whether to delete the compressed datasets.
#' If FALSE, the tar files will be stored in a folder called datasets_tar
#' under the \code{output_dir}.
#' @param quiet logical whether to display no messages.
#'
#' @author Dongze He
#'
#'
#' @details
#' 10x Genomics provided various publicly available datasets on
#' their website for free downloading
#' (\url{https://www.10xgenomics.com/resources/datasets}).
#' To stop re-inventing the wheel, we downloaded and processed
#' these datasets using a Nextflow-based alevin-fry workflow
#' (\url{https://github.com/COMBINE-lab/10x-requant}) and
#' provide the link to the quantification results in the above
#' GitHub repository. Using this function, one can directly access
#' the quantification results of the preprocessed 10x datasets
#' and load these results into R as a SingleCellExperiment object.
#' Currently, the datasets that are available for querying include:
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
#' To obtain these information as a dataframe, one can simply run
#' `preprocessed_10x_data()` in R.
#' Note that because the name of datasets are too long, the stored
#' datasets are named by the MD5 hash value of their fastqs.tar file. 
#' If one would like to use the downloaded quantification results
#' out of this funbction, please refer to the webpage of the
#' corresponding datasets to find the MD5 hash value.
#'
#' @export
#'
#' @return If return_available_dataset_df is set as true,
#' a dataframe containing all available datasets will be returned;
#' otherwise, a vector of local paths to the downloaded datasets
#' will be returned.
#'
#' @examples
#' 
#' \dontrun{
#' library(roe)
#' # run the function
#' available_datasets = preprocessed_10x_data()
#' preprocessed_10x_data(dataset_id = c(1, 2),
#'                        output_dir = "10x_datasets",
#'                        force = FALSE,
#'                        delete_tar = TRUE,
#'                        quiet = FALSE
#' )
#' }
#'



preprocessed_10x_data <- function(dataset_ids = c(),
                              output_dir = "10x_datasets",
                              force = FALSE,
                              delete_tar = TRUE,
                              quiet = FALSE
                              ) {

# available_datasets = read.csv("available_datasets.tsv", sep = "\t")
# usethis::use_data(available_datasets, internal = TRUE)

  # if the user just wants the data frame, return it
  if (length(dataset_ids) == 0) {
    return (available_datasets)
  }

  .say(quiet, "Check the validity of dataset_ids")

  # now check the validity of dataset_ids
  for (idx in seq(length(dataset_ids))) {
    dataset_id <- dataset_ids[idx]
    if (!(is.numeric(dataset_id) & dataset_id <= nrow(available_datasets))) {
      message("Found invalid dataset id: ", dataset_id, ", ignored")
      dataset_ids <- dataset_ids[-idx]
    }
  }

  # check whether there is any dataset id left
  if (length(dataset_ids) == 0) {
    stop("No valid dataset id found, can not proceed")
  }

  # download the quantification tar file for each queried dataset.
  quant_dir_list <- c()
  # folder for (temporarily) storing tar files.
  tar_dir <- file.path(output_dir, "datasets_tar")
  dir.create(tar_dir, recursive = TRUE,
              showWarnings = FALSE)

  for (dataset_id in dataset_ids) {
    .say(quiet, "\n\nProceeding dataset #", dataset_id)

    # specify paths
    quant_parent_dir <- file.path(output_dir,
                          available_datasets[dataset_id, "MD5"]
                        )

    tar_file <- file.path(tar_dir,
                            paste0(available_datasets[dataset_id, "MD5"],
                                  ".tar")
                                  )

    # process if needed
    if (file.exists(quant_parent_dir)) {
      .say(quiet, "    output dir exists: \n    ", quant_parent_dir, "\n")

      if (force) {
        .say(quiet, "    force re-processing")
        unlink(quant_parent_dir, recursive = TRUE, force = TRUE)
      } else {
        quant_dir <- list.dirs(quant_parent_dir,
                                full.names = TRUE,
                                recursive = FALSE)
        quant_dir_list <- c(quant_dir_list, quant_dir)
        next
      }
    }

    .say(quiet, "    Downloading alevin-fry quant folder")
    url <- available_datasets[dataset_id, "quant_link"]
    utils::download.file(url = url,
                          destfile = tar_file,
                          quiet = TRUE,
                          cacheOK = FALSE)

    .say(quiet, "    Decompressing alevin-fry quant folder")
    utils::untar(tarfile = tar_file,
                  exdir = quant_parent_dir
                )
    quant_dir <- list.dirs(quant_parent_dir,
                            full.names = TRUE,
                            recursive = FALSE)

    quant_dir_list <- c(quant_dir_list, quant_dir)
    .say(quiet, "\n")

  }

  if (delete_tar) {
    .say(quiet, "Delete temp tar files")
    unlink(tar_dir,  recursive = TRUE, force = TRUE)
  }
  if (!quiet) {
    message("Done")
  }
  names(quant_dir_list) <- dataset_ids
  quant_dir_list
}
