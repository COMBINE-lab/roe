#' @export
#' @rdname ProcessedQuant
#' 
methods::setClass("ProcessedQuant",
            methods::representation(dataset_id = "numeric",
                                    chemistry = "character",
                                    reference = "character",
                                    dataset_name = "character",
                                    dataset_url = "character",
                                    fastq_url = "character",
                                    fastq_MD5sum = "character",
                                    feature_barcode_csv_url = "character",
                                    multiplexing_library_csv_url = "character",
                                    quant_tar_url = "character",
                                    quant_path = "character",
                                    tar_path = "character",
                                    sce = "SingleCellExperiment")
                )

#' return the dataframe that contains the details of all
#' available datasets.
#' @export

get_available_dataset_df <- function() {
    # available_datasets = read.csv("available_datasets.tsv", sep = "\t")
    # usethis::use_data(available_datasets, internal = TRUE, force = TRUE)
    available_datasets
}

#' print the index and the name of available datasets.
#'
#' @export
print_available_datasets <- function() {
    for (dataset_id in seq(nrow(available_datasets))) {
        cat(paste0(dataset_id,
                    ". ",
                    available_datasets[dataset_id, "dataset_name"],
                    "\n"))
    }
}

#' The ProcessedQuant class
#' The ProcessedQuant class is designed to represent the details of
#' preprocessed dataset from the
#' \href{https://github.com/COMBINE-lab/10x-requant}{10x-requant}
#' project. To initialize a ProcessedQuant class, run
#' \code{ProcessedQuant(dataset_id)} where \code{dataset_id} is
#' the id of an available dataset. Run \code{print_available_datasets()}
#' for obtaining the id of all available datasets.
#' The instantiated ProcessedQuant object needs to go through some
#' function calls, i.e, \code{fetch_quant()}, \code{decompress_quant()}
#' and \code{load_quant()} to become a complete ProcessedQuant class,
#' in which the \code{sce} slot contains the
#' \code{SingleCellExperiment} object of the quantification result
#' of the dataset. information of a specific dataset.
#' @param dataset_id the id of an available dataset.
#' 
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment

ProcessedQuant <- function(dataset_id) {
    if ((is.numeric(dataset_id))) {
        if (!((dataset_id >= 1) & (dataset_id <= nrow(available_datasets)))) {
            stop("Invalid dataset_id, run print_available_datasets()",
                " to get available dataset ids.")
        }
    } else {
        stop("Invalid dataset_id type, accepts integer only.")
    }
    dataset_info <- available_datasets[dataset_id, ]
    processed_quant <- methods::new(
                            "ProcessedQuant",
                            dataset_id = dataset_info$"dataset_id",
                            chemistry = dataset_info$"chemistry",
                            reference = dataset_info$"reference",
                            dataset_name = dataset_info$"dataset_name",
                            dataset_url = dataset_info$"dataset_url",
                            fastq_url = dataset_info$"fastq_url",
                            fastq_MD5sum = dataset_info$"fastq_MD5sum",
                            quant_tar_url = dataset_info$"quant_tar_url",
                            sce = SingleCellExperiment()
    )

    if (!is.na(dataset_info$feature_barcode_csv_url)) {
        processed_quant@feature_barcode_csv_url <-
            dataset_info$feature_barcode_csv_url
    }
    if (!is.na(dataset_info$multiplexing_library_csv_url)) {
        processed_quant@multiplexing_library_csv_url <-
            dataset_info$multiplexing_library_csv_url
    }

    processed_quant
}

#' fetch the compressed quantification result of a
#' specific dataset and stroe the path to the tar file
#' in the tar_path slot of the returned processed_quant list.
#' This function must be run after \code{ProcessedQuant()}
#' 
#' @param processed_quant a \code{ProcessedQuant} class object
#' recording the details of a dataset, which is nitialized by
#' running \code{ProcessedQuant(dataset_id)}
#' @param tar_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param tar_name a string indicates the name of the tar file. If NULL,
#' the dataset_id will be used as the file name.
#' @param force logic whether to proceed if the tar file exists.
#' @param quiet logical whether to display no messages
#' @return A ProcessedQuant class object with a updated \code{tar_path} slot.
#' @export

fetch_quant <- function(processed_quant,
                        tar_dir="quant_tar",
                        tar_name=NULL,
                        force=FALSE,
                        quiet=FALSE) {
    # check validity of processed_quant
    check_validity(processed_quant)
    .say(quiet,
        "Fetching the quant result of dataset #",
        processed_quant@dataset_id
    )

    # if tar_path is not null, return unless force=TRUE
    if (!identical(processed_quant@tar_path, character(0))) {
        if (file.exists(processed_quant@tar_path)) {
            if (!force) {
                .say(quiet,
                    "  - The processed_quant@tar_path is not ",
                    "empty and the path exists:\n",
                    "    ", processed_quant@tar_path, "\n",
                    "  - Pass force=TRUE to update it")
                return(processed_quant)
            }
        }
    }

    dir.create(tar_dir, recursive = TRUE,
                showWarnings = FALSE)

    if (is.null(tar_name)) {
        tar_name <- paste0(processed_quant@dataset_id, ".tar")
    } else if (!endsWith(tar_name, ".tar")) {
        tar_name <- paste0(tar_name, ".tar")
    }

    tar_file <- file.path(tar_dir, tar_name)

    if (file.exists(tar_file)) {
        if (force) {
            .say(quiet,
                "  - Overwriting the existing tar file:\n",
                "    ", tar_file, "\n")
        } else {
            .say(quiet,
                "  - Use the existing file as tar_path:\n",
                "    ", tar_file, "\n",
                "  - Pass force=TRUE to overwrite it")
            processed_quant@tar_path = tar_file
            return(processed_quant)
        }
    }

    url <- processed_quant@quant_tar_url
    utils::download.file(url = url,
                        destfile = tar_file,
                        quiet = TRUE,
                        cacheOK = FALSE
    )
    processed_quant@tar_path <- tar_file

    .say(quiet,
        "  - Fetched quant tar is saved as:\n",
        "    ", processed_quant@tar_path
    )
    return(processed_quant)
}

#' decompress the fetched quantification result of a
#' specific dataset and record the path
#' as the quant_path slot of the returned processed_quant list
#' This function must be run after \code{fetch_quant()}
#' 
#' @param processed_quant a \code{ProcessedQuant} class object
#' recording the details of
#' a dataset. Initialized by \code{ProcessedQuant(dataset_id)}
#' @param quant_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param quant_path_name a string indicates the name of the directory
#' that will be used for storing the quantification result. If NULL,
#' the dataset_id will be used as the file name.
#' @param force logic whether to proceed if the tar file exists.
#' @param quiet logical whether to display no messages
#' @return A ProcessedQuant class object with a updated \code{quant_path} slot.
#' @export

decompress_quant <- function(processed_quant,
                            quant_dir="processed_quant",
                            quant_path_name=NULL,
                            force=FALSE,
                            quiet=FALSE) {
    check_validity(processed_quant)

    if (identical(processed_quant@tar_path, character(0))) {
        stop("tar_path slot is empty, ",
            "run processed_quant = fetch_quant(processed_quant) ",
            "to fetch the tar file.")
    }

    .say(quiet,
        "Decompressing the quant result of dataset #",
        processed_quant@dataset_id, " using: \n",
        "  ", processed_quant@tar_path
    )

    # if quant_path is not null, return unless force=TRUE
    if (!identical(processed_quant@quant_path, character(0))) {
        if (file.exists(processed_quant@quant_path)) {
            if (!force) {
                .say(quiet,
                    "  - Use the existing directory as quant_path:\n",
                    "    ", processed_quant@quant_path, "\n",
                    "  - Pass force=TRUE to update it")
                return(processed_quant)
            }
        }
    }

    # check quant_path_name
    if (is.null(quant_path_name)) {
        quant_path_name <- paste0(processed_quant@dataset_id)
    }

    # specify paths
    quant_parent_dir <- file.path(quant_dir, quant_path_name)

    # check quant_parent_dir
    if (file.exists(quant_parent_dir)) {
        if (force) {
            .say(quiet,
                "  - Removing existing quant folder:\n",
                "    ", quant_parent_dir)
            unlink(quant_parent_dir, recursive = TRUE, force = TRUE)
        } else {
            processed_quant@quant_path <- list.dirs(quant_parent_dir,
                                                    full.names = TRUE,
                                                    recursive = FALSE)
            .say(quiet,
                "  - Use the existing directory as quant_path:",
                "    ", processed_quant@quant_path,
                "  - pass force=TRUE to overwrite it")
            return(processed_quant)
        }
    }

    # if we are here, untar it
    utils::untar(tarfile = processed_quant@tar_path,
                exdir = quant_parent_dir
    )
    processed_quant@quant_path <- list.dirs(quant_parent_dir,
                                            full.names = TRUE,
                                            recursive = FALSE)

    .say(quiet,
        "  - Decompressed quant result is saved as:\n",
        "    ", processed_quant@quant_path
    )
    return(processed_quant)
}

#' load the fetched quantification result of a
#' specific dataset as a SingleCellExperiment object
#' and store it in the sce slot of the returned
#'  processed_quant list as the sce slot.
#' This function must be run after \code{decompress_quant()}
#' 
#' @param processed_quant a \code{ProcessedQuant} class object  
#' recording the details of
#' a dataset. 
#' @param output_format the format of the returned SingleCellExperiment
#' object. It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{outputFormat} parameter.
#' @param nonzero It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{nonzero} parameter.
#' @param force logic whether to proceed if the sce slot exists.
#' @param quiet logical whether to display no messages. Default is set as FALSE.
#' @return A ProcessedQuant class object with a updated \code{sce} slot.
#' 
#' @export

load_quant <- function(processed_quant,
                        output_format="scRNA",
                        nonzero = FALSE,
                        force = TRUE,
                        quiet = FALSE) {
    check_validity(processed_quant)

    if (identical(processed_quant@quant_path, character(0))) {
        stop("quant_path slot is invalid, ",
            "run processed_quant= dec",
            "ompress_tar(processed_quant)",
            "to prepare it.")
    }

    if (!file.exists(processed_quant@quant_path)) {
        stop("quant_path slot is invalid, ",
                "run processed_quant= dec",
                "ompress_tar(processed_quant)",
                "to prepare it.")
    }

    if ((sum(dim(processed_quant@sce)) != 0) & (!force)) {
        .say(quiet,
            "  - The sce slot of the passed processed_",
            "quant list is not empty",
            "  - Pass force=TRUE to update it\n")
        return(processed_quant)
    }

    .say(quiet,
        "Loading dataset #",
        processed_quant@dataset_id,
        " from:\n  ",
        processed_quant@quant_path
    )
    processed_quant@sce <- fishpond::loadFry(
                                fryDir = processed_quant@quant_path,
                                outputFormat = output_format,
                                nonzero = nonzero,
                                quiet = quiet
                            )

    return(processed_quant)
}

#' fetch, decompress and load the quantification result of a
#' specific dataset as a SingleCellExperiment object
#' and store it in the sce slot of the returned
#' processed_quant list.
#' 
#' @param dataset_id the id of an available dataset.
#' Run \code{print_available_datasets()} for all available
#' datasets.
#' @param tar_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param tar_name a string indicates the name of the tar file. If NULL,
#' the dataset_id will be used as the file name.
#' @param quant_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param quant_path_name a string indicates the name of the directory
#' that will be used for storing the quantification result. If NULL,
#' the dataset_id will be used as the file name.
#' @param output_format the format of the returned SingleCellExperiement
#' object. It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{outputFormat} parameter.
#' @param force logic whether to proceed if the files exist.
#' @param nonzero It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{nonzero} parameter.
#' @param quiet logical whether to display no messages
#' @return A complete ProcessedQuant class object with valid
#' \code{tar_path}, \code{quant_path} and \code{sce} slots.
#' @export
#' 
FDL <- function(dataset_id,
                tar_dir="quant_tar",
                tar_name=NULL,
                quant_dir="processed_quant",
                quant_path_name=NULL,
                output_format="scRNA",
                nonzero=FALSE,
                force=FALSE,
                quiet=FALSE) {

    # init processed_quant
    processed_quant <- ProcessedQuant(dataset_id)

    # fetch it
    processed_quant <- fetch_quant(processed_quant,
                                    tar_dir = tar_dir,
                                    tar_name = tar_name,
                                    force = force,
                                    quiet = quiet)

    # decompress it
    processed_quant <- decompress_quant(processed_quant,
                                        quant_dir = quant_dir,
                                        quant_path_name = quant_path_name,
                                        force = force,
                                        quiet = quiet)

    # load it
    processed_quant <- load_quant(processed_quant,
                                    output_format = output_format,
                                    nonzero = nonzero,
                                    quiet = quiet)

    return(processed_quant)
}

check_validity <- function(processed_quant) {
    if (any(c(
        identical(processed_quant@dataset_id, numeric(0)),
        identical(processed_quant@chemistry, character(0)),
        identical(processed_quant@reference, character(0)),
        identical(processed_quant@dataset_name, character(0)),
        identical(processed_quant@dataset_url, character(0)),
        identical(processed_quant@fastq_url, character(0)),
        identical(processed_quant@fastq_MD5sum, character(0)),
        identical(processed_quant@quant_tar_url, character(0))
    ))) {
        stop(paste0("Invalid processed_quant, use ",
                    "ProcessedQuant(dataset_id) to ",
                    "initiate a valid one.")
        )
    } else {
        dataset_id <- processed_quant@dataset_id
        if ((is.numeric(dataset_id))) {
            if (!((dataset_id >= 1) &
                (dataset_id <= nrow(available_datasets)))) {
                stop(paste0("Invalid dataset info list, use ",
                            "ProcessedQuant(dataset_id) to ",
                            "initiate a valid one.")
                )
            }
        } else {
            stop(paste0("Invalid dataset info list, use ",
                        "ProcessedQuant(dataset_id) to ",
                        "initiate a valid one.")
            )
        }
    }
}
