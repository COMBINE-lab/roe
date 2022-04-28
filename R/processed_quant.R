get_available_dataset_df <- function() {
    # available_datasets = read.csv("available_datasets.tsv", sep = "\t")
    # usethis::use_data(available_datasets, internal = TRUE, force = TRUE)
    available_datasets
}

print_available_datasets <- function() {
    for (dataset_id in seq(nrow(available_datasets))) {
        cat(paste0(dataset_id,
                    ". ",
                    available_datasets[dataset_id, "dataset_name"],
                    "\n"))
    }
}

init_dataset_info_list <- function(dataset_id) {
    if ((is.numeric(dataset_id))) {
        if (!((dataset_id >= 1) & (dataset_id <= nrow(available_datasets)))) {
            stop(paste0("Invalid dataset_id, run print_available_datasets()",
                        " to get available dataset ids."))
        }
    } else {
        stop(paste0("Invalid dataset_id type, accepts integer only."))
    }
    dataset_info <- available_datasets[dataset_id, ]
    dataset_info_list <- list(dataset_id = dataset_info$"dataset_id",
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
    dataset_info_list
}


fetch_tar <- function(dataset_info_list,
                        tar_dir="quant_tar",
                        file_name=NULL,
                        force=FALSE,
                        quiet=FALSE) {
    # check validity of dataset_info_list
    check_validity(dataset_info_list)
    .say(quiet, paste0("Fetching the quant result of dataset #",
                        dataset_info_list[["dataset_id"]]))

    if ((!is.null(dataset_info_list[["tar_path"]])) &
        file.exists(dataset_info_list[["tar_path"]])) {
        .say(quiet, "  - The tar_path field is not None and the path exists:")
        .say(quiet, "    ", dataset_info_list[["tar_path"]], "\n")
        .say(quiet, "  - Pass force=True to fetch it again\n")
        return(dataset_info_list)
    }

    dir.create(tar_dir, recursive = TRUE,
                showWarnings = FALSE)

    if (is.null(file_name)) {
        file_name <- paste0(dataset_info_list[["dataset_id"]], ".tar")
    } else if (!endsWith(file_name, ".tar")) {
        file_name <- paste0(file_name, ".tar")
    }

    tar_file <- file.path(tar_dir, file_name)

    .say(quiet, "    Downloading alevin-fry quant folder")
    url <- dataset_info_list[["quant_link"]]
    utils::download.file(url = url,
                        destfile = tar_file,
                        quiet = TRUE,
                        cacheOK = FALSE)

    .say(quiet, "    Decompressing alevin-fry quant folder")
    utils::untar(tarfile = tar_file,
                exdir = quant_parent_dir
                )
    dataset_info_list[["quant_dir"]] <- list.dirs(quant_parent_dir,
                                                full.names = TRUE,
                                                recursive = FALSE)

    dataset_info_list_list[[as.character(dataset_id)]] <- dataset_info_list_list
    .say(quiet, "\n")
}

check_validity <- function(dataset_info_list) {
    if (any(c(
        is.null(dataset_info_list[["dataset_id"]]),
        is.null(dataset_info_list[["chemistry"]]),
        is.null(dataset_info_list[["reference"]]),
        is.null(dataset_info_list[["dataset_name"]]),
        is.null(dataset_info_list[["link"]]),
        is.null(dataset_info_list[["data_url"]]),
        is.null(dataset_info_list[["MD5"]]),
        is.null(dataset_info_list[["feature_barcode"]]),
        is.null(dataset_info_list[["library_csv"]]),
        is.null(dataset_info_list[["quant_link"]])
    ))) {
        stop(paste0("Invalid dataset info list, use init_",
                    "dataset_info_list(dataset_id) to ",
                    "initiate a valid one.")
            )
    } else {
        dataset_id = dataset_info_list[["dataset_id"]]
        if ((is.numeric(dataset_id))) {
            if (!((dataset_id >= 1) &
                    (dataset_id <= nrow(available_datasets)))) {
                        stop(paste0("Invalid dataset info list, use init_",
                                    "dataset_info_list(dataset_id) to ",
                                    "initiate a valid one.")
                                    )
            }
        } else {
            stop(paste0("Invalid dataset info list, use init_",
                        "dataset_info_list(dataset_id) to ",
                        "initiate a valid one.")
                        )
        }
    }
}
