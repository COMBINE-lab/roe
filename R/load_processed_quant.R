#' Fetch and load 10x datasets
#'
#' Fetch and load preprocessed 10x datasets as SingleCellExperiment objects.
#' @param dataset_ids integer scalar or vector providing the id of the
#' dataset(s) to be Fetched and processed. See \code{\link{fetch_processed_quant}}
#' for details.
#' @param fetch_dir path to the folder where the fetched quantification results
#' should be stored.
#' See \code{\link{fetch_processed_quant}} for details.
#' @param force logical whether to force re-downloading the existing datasets.
#' See \code{\link{fetch_processed_quant}} for details.
#' @param delete_tar logical whether to delete the intermediate 
#' compressed datasets.
#' See \code{\link{fetch_processed_quant}} for details.
#' @param output_format can be \emph{either} a valid
#' \code{\link[fishpond]{loadFry}} \code{outputFormat}
#' or a list of it. If it is a list, each item in the list
#' corresponds to a queried dataset and
#' should be a valid loadFry-\code{outputFormat}.
#' Each queried dataset will be loaded according to
#' the corresponding \code{outputFormat} in the list.
#' @param nonzero Similar with \code{outputFormat}, it can be
#' either a valid \code{\link[fishpond]{loadFry}} \code{nonzero}
#' parameter or a list of it. If a list, the queried datasets
#' will be processed using the corresponding \code{nonzero}.
#' @param quiet logical whether to display no messages
#'
#' @details
#' This function is a thin wrapper of the \code{\link[fishpond]{loadFry}}
#' function and \code{\link{fetch_processed_quant}} function. This function will
#' fetch the quantification result of each dataset in the provided 
#' `dataset_ids` into the `fetch_dir` directory, and load the 
#' quantification result of each fetched dataset into R as a 
#' SingleCellExperiment object. This function returns a list of lists, 
#' in which each list stores the information of a fetched dataset, and 
#' the SingleCellExperiment of that dataset.
#' To list all available datasets,
#' simply run \code{load_processed_quant()}.
#' This function takes the combination
#' of the parameters of the \code{\link[fishpond]{loadFry}}
#' function and \code{\link{fetch_processed_quant}} function. For
#' each of the \code{\link[fishpond]{loadFry}} parameters,
#' a list of that parameter can be specified, each for a
#' fetched dataset.
#' 
#' @export
#'
#' @return If an empty dataset_ids is provided,
#' a data frame containing the information of available datasets will be returned;
#' otherwise, a list of lists, where each list stores the information of one 
#' fetched dataset. The `quant_dir` field represent the path to the quantification
#' result of the fetched dataset; The `sce` field represent the SingleCellExperiment
#' of this dataset.
#'
#' @examples
#' 
#' \dontrun{
#' library(roe)
#' # run the function
#' # The four way to define output_format
#' # are the same.
#' load_processed_quant(dataset_ids = c(1, 2),
#'         fetch_dir = "10x_datasets",
#'         force = FALSE,
#'         delete_tar = TRUE,
#'         output_format = "scRNA",
#' #         output_format = list("scRNA", "scRNA"),
#' #         output_format = list(list(counts = c("S", "A")),
#' #                               list(counts = c("S", "A"))
#' #                              ),
#' #         output_format = list(counts = c("S", "A")),
#'         nonzero = FALSE,
#'         quiet = FALSE
#' )
#' }
#' 

load_processed_quant <- function(dataset_ids = c(),
                    fetch_dir = "processed_quant",
                    force = FALSE,
                    delete_tar = TRUE,
                    output_format = "scRNA",
                    nonzero = FALSE,
                    quiet = FALSE
) {
    .say(quiet, "Processing parameters")
    # check whether output_format are valid
    # we just check the length, the validity of
    # each outputFormat will be checked by loadFry
    nd <- length(dataset_ids)
    
    # if the user just wants the data frame, return it
    if (length(dataset_ids) == 0) {
        return (available_datasets)
    }

    if (is.list(output_format)) {
        # if a list is given,
        # it should be either one customized format
        # or the format of each fetched datasets
        # so check the name
        if (!setequal(sort(names(output_format)), sort(dataset_ids))) {
            # now it should be one customized format
            output_format <- as.list(rep(list(output_format), nd))
            names(output_format) <- dataset_ids
        }
        # otherwise, each dataset should get a format, so check the length
        if (length(output_format) != nd) {
            stop("The providing output_format list has different length with dataset_ids, cannot proceed")
        }
    } else if (is.character(output_format)) {
        # if a str is given, it should be a pre-defined format
        # and it will be used for all datasets
        output_format <- as.list(rep(output_format, nd))
        names(output_format) <- dataset_ids
    } else {
        stop("The providing output_format list has different length with dataset_ids, cannot proceed")
    }

    # then we do the same thing for nonzero
    if (is.list(nonzero)) {
        # user defines a nonzero for each dataset
        # check whether length matches
        if (length(nonzero) != nd) {
            stop("The providing nonzero list has different length with dataset_ids, cannot proceed")
        }
    } else {
        # The last valid situation is that the
        # user provides a sinlg pre-defined format
        # then we repeat it for every dataset
        nonzero <- as.list(rep(nonzero, nd))
        names(nonzero) <- dataset_ids
    }

    .say(quiet, "Fetching datasets")

    # download the datsets
    dataset_paths <- fetch_processed_quant(dataset_ids = dataset_ids,
                                            fetch_dir = fetch_dir,
                                            force = force,
                                            delete_tar =  delete_tar,
                                            quiet = quiet
                                            )

    processed_dataset_list <- list()
    # process them using user output
    for (dataset_id in dataset_ids) {
        dataset_id <- as.character(dataset_id)
        .say(quiet, "Loading dataset ", dataset_id)
        dataset_path_ds <- dataset_paths[dataset_id]
        output_format_ds <- output_format[[dataset_id]]
        nonzero_ds <- nonzero[[dataset_id]]
        
        processed_dataset = .get_dataset_info_list(available_datasets, dataset_id)
        processed_dataset[["quant_dir"]] = dataset_paths[dataset_id]
        processed_dataset[[dataset_id]] <-
                            fishpond::loadFry(fryDir = dataset_path_ds,
                                                outputFormat = output_format_ds,
                                                nonzero = nonzero_ds,
                                                quiet = quiet
                        )
        processed_dataset_list[[dataset_id]] = processed_dataset
    }

    # output
    processed_dataset_list
}
