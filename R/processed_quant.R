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


#' initiate the processed_quant list used for recording the 
#' information of a specific dataset
#' @param dataset_id the id of an available dataset. 
#' Run \code{print_available_datasets()} for all available
#' datasets.
#' 
#' @export

init_processed_quant <- function(dataset_id) {
    if ((is.numeric(dataset_id))) {
        if (!((dataset_id >= 1) & (dataset_id <= nrow(available_datasets)))) {
            stop("Invalid dataset_id, run print_available_datasets()",
                " to get available dataset ids.")
        }
    } else {
        stop("Invalid dataset_id type, accepts integer only.")
    }
    dataset_info <- available_datasets[dataset_id, ]
    processed_quant <- list(dataset_id = dataset_info$"dataset_id",
                            chemistry = dataset_info$"chemistry",
                            reference = dataset_info$"reference",
                            dataset_name = dataset_info$"dataset_name",
                            link = dataset_info$"link",
                            data_url = dataset_info$"data_url",
                            MD5	= dataset_info$"MD5",
                            feature_barcode	= dataset_info$"feature_barcode",
                            library_csv	= dataset_info$"library_csv",
                            quant_link = dataset_info$"quant_link",
                            quant_path = NULL,
                            tar_path = NULL,
                            sce = NULL
    )
    processed_quant
}

#' fetch the compressed quantification result of a 
#' specific dataset and stroe the path to the tar file 
#' in the tar_path field of the returned processed_quant list.
#' This function must be run after \code{init_processed_quant()}
#' 
#' @param processed_quant a list recording the details of 
#' a dataset. Initialized by \code{init_processed_quant(dataset_id)}
#' @param tar_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param file_name a string indicates the name of the tar file. If NULL,
#' the dataset_id will be used as the file name.
#' @param force logic whether to proceed if the tar file exists.
#' @param quiet logical whether to display no messages
#' @export

fetch_tar <- function(processed_quant,
                        tar_dir="quant_tar",
                        file_name=NULL,
                        force=FALSE,
                        quiet=FALSE) {
    # check validity of processed_quant
    check_validity(processed_quant)
    .say(quiet, 
         "Fetching the quant result of dataset #",
         processed_quant$dataset_id
        )

    # if tar_path is not null, return unless force=TRUE
    if (file.exists(processed_quant$tar_path)) {
        .say(quiet,
             "  - The processed_quant$tar_path exists:\n",
             "    ", processed_quant$tar_path, "\n",
             "  - Pass force=True to update it\n")
        return(processed_quant)
    }

    dir.create(tar_dir, recursive = TRUE,
                showWarnings = FALSE)

    if (is.null(file_name)) {
        file_name <- paste0(processed_quant$dataset_id, ".tar")
    } else if (!endsWith(file_name, ".tar")) {
        file_name <- paste0(file_name, ".tar")
    }    
    
    tar_file <- file.path(tar_dir, file_name)

    if (file.exists(tar_file)) {
        if (force) {
            .say(quiet,
                 "  - Overwriting the existing tar file:\n",
                 "    ", tar_file, "\n")
        } else {
            .say(quiet,
                 "  - Use the existing file as tar_path:\n",
                 "    ", tar_file, "\n",
                 "  - Pass force=True to overwrite it\n")
            processed_quant$tar_path = tar_file
            return(processed_quant)
        }
    }

    url <- processed_quant$quant_link
    utils::download.file(url = url,
                        destfile = tar_file,
                        quiet = TRUE,
                        cacheOK = FALSE
                        )
    processed_quant$tar_path = tar_file
    
    say(quiet, 
        "  - Fetched quant tar is saved as:\n",
        "    ", processed_quant$tar_path
        )
}

#' decompress the fetched quantification result of a 
#' specific dataset and record the path 
#' as the quant_path field of the returned processed_quant list
#' This function must be run after \code{fetch_tar()}
#' 
#' @param processed_quant a list recording the details of 
#' a dataset. Initialized by \code{init_processed_quant(dataset_id)}
#' @param quant_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param quant_path_name a string indicates the name of the directory 
#' that will be used for storing the quantification result. If NULL,
#' the dataset_id will be used as the file name.
#' @param force logic whether to proceed if the tar file exists.
#' @param quiet logical whether to display no messages
#' @export

decompress_tar <- function(processed_quant,
                           quant_dir="processed_quant",
                           quant_path_name=None,
                           force=FALSE,
                           quiet=FALSE) {
    check_validity(processed_quant)
    
    if (is.null(processed_quant$tar_path)) {
        stop("tar_path field is NULL, ",
             "run processed_quant = fetch_tar(processed_quant) ",
             "to fetch the tar file.")
    }
    
    .say(quiet, 
         "Decompressing the quant result of dataset #",
         processed_quant$dataset_id
    )
    
    # if quant_path is not null, return unless force=TRUE
    if (file.exists(processed_quant$quant_path)) {
        .say(quiet,
             "  - Use the existing directory as quant_path:\n",
             "    ", processed_quant$quant_path, "\n",
             "  - Pass force=True to update it\n")
        return(processed_quant)
    }
    
    # check quant_path_name
    if (is.null(quant_path_name)) {
        quant_path_name = paste0(processed_quant$dataset_id)
    }
    
    # specify paths
    quant_parent_dir <- file.path(fetch_dir, quant_path_name)

    # check quant_parent_dir
    if (file.exists(quant_parent_dir)) {
        if (force) {
            .say(quiet,
                 "  - Removing existing quant folder:\n",
                 "    ", quant_parent_dir)
            unlink(quant_parent_dir, recursive = TRUE, force = TRUE)
        } else {
            processed_quant$quant_path <- list.dirs(quant_parent_dir,
                                                           full.names = TRUE,
                                                           recursive = FALSE)
            say(quiet, 
                "  - Use the existing directory as quant_path:",
                "    ", processed_quant$quant_path,
                "  - pass force=True to overwrite it\n")
            return(processed_quant)
        }
    }
    
    # if we are here, untar it
    utils::untar(tarfile = tar_file,
                 exdir = quant_parent_dir
    )
    processed_quant$quant_path <- list.dirs(quant_parent_dir,
                                                  full.names = TRUE,
                                                  recursive = FALSE)
    
    return(processed_quant)
}

#' load the fetched quantification result of a 
#' specific dataset as a SingleCellExperiment object
#' and store it in the sce field of the returned
#'  processed_quant list.
#' This function must be run after \code{decompress_tar()}
#' 
#' @param processed_quant a list recording the details of 
#' a dataset. 
#' @param output_format the format of the returned SingleCellExperiement
#' object. It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{outputFormat} parameter. 
#' @param nonzero It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{nonzero} parameter. 
#' @param quiet logical whether to display no messages
#' @export

load_quant <- function(processed_quant,
                       output_format="scRNA",
                       nonzero = FALSE,
                       quiet = FALSE) {
    check_validity(processed_quant)
    
    if (!file.exists(processed_quant$quant_path)) {
        stop("quant_path field is invalid, ",
             "run processed_quant= dec",
             "ompress_tar(processed_quant)",
             "to prepare it.")
    }

    .say(quiet, 
         "Decompressing the quant result of dataset #",
         processed_quant$dataset_id,
         " from: ",
         processed_quant$quant_path
    )
    processed_quant$sce <- fishpond::loadFry(fryDir = processed_quant$quant_path,
                      outputFormat = output_format,
                      nonzero = nonzero,
                      quiet = quiet)
    
    return(processed_quant)
}

#' fetch, decompress and load the quantification result of a 
#' specific dataset as a SingleCellExperiment object
#' and store it in the sce field of the returned
#' processed_quant list.
#' 
#' @param dataset_id the id of an available dataset. 
#' Run \code{print_available_datasets()} for all available
#' datasets.
#' @param tar_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param tar_file_name a string indicates the name of the tar file. If NULL,
#' the dataset_id will be used as the file name.
#' @param quant_dir a string to a path where the fetched tar files will
#' be stored. It will be created if does not exist.
#' @param quant_path_name a string indicates the name of the directory 
#' that will be used for storing the quantification result. If NULL,
#' the dataset_id will be used as the file name.
#' @param output_format the format of the returned SingleCellExperiement
#' object. It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{outputFormat} parameter. 
#' @param nonzero It will be passed to \code{\link[fishpond]{loadFry}}
#' as the \code{nonzero} parameter. 
#' @param quiet logical whether to display no messages
#' @export
#' 
FDL <- function(dataset_id,
                tar_dir="quant_tar",
                tar_file_name=NULL,
                quant_dir="processed_quant",
                quant_path_name=NULL,
                output_format="scRNA",
                nonzero=FALSE,
                force=FALSE, 
                quiet=FALSE) {

    # init processed_quant
    processed_quant = init_processed_quant(dataset_id)
    
    # fetch it
    processed_quant = fetch_tar(processed_quant,
                              tar_dir=tar_dir,
                              file_name=tar_file_name,
                              force=force,
                              quiet=quiet)
    
    # decompress it
    processed_quant = decompress_tar(processed_quant,
                                   quant_dir=quant_dir,
                                   quant_path_name=quant_path_name,
                                   force=force,
                                   quiet=quiet)
    
    # load it
    processed_quant = load_quant(processed_quant,
                               output_format=output_format,
                               nonzero = nonzero,
                               quiet = quiet)
    
    return(processed_quant)
}

check_validity <- function(processed_quant) {
    if (any(c(
        is.null(processed_quant$dataset_id),
        is.null(processed_quant$chemistry),
        is.null(processed_quant$reference),
        is.null(processed_quant$dataset_name),
        is.null(processed_quant$link),
        is.null(processed_quant$data_url),
        is.null(processed_quant$MD5),
        is.null(processed_quant$feature_barcode),
        is.null(processed_quant$library_csv),
        is.null(processed_quant$quant_link)
    ))) {
        stop(paste0("Invalid dataset info list, use init_",
                    "processed_quant(dataset_id) to ",
                    "initiate a valid one.")
            )
    } else {
        dataset_id = processed_quant$dataset_id
        if ((is.numeric(dataset_id))) {
            if (!((dataset_id >= 1) &
                    (dataset_id <= nrow(available_datasets)))) {
                        stop(paste0("Invalid dataset info list, use init_",
                                    "processed_quant(dataset_id) to ",
                                    "initiate a valid one.")
                                    )
            }
        } else {
            stop(paste0("Invalid dataset info list, use init_",
                        "processed_quant(dataset_id) to ",
                        "initiate a valid one.")
                        )
        }
    }
}
