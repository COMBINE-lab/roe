#' Download and load 10x datasets
#'
#' Download and load preprocessed 10x datasets as SingleCellExperiment objects.
#' @param dataset_ids integer scalar or vector providing the id of the
#' dataset(s) to be downloaded and processed. See \code{\link{preprocessed_10x_data}}
#' for details.
#' @param output_dir path to the folder that will
#' be used to store downloaded datasets.
#' See \code{\link{preprocessed_10x_data}} for details.
#' @param force logical whether to force re-downloading the existing datasets.
#' See \code{\link{preprocessed_10x_data}} for details.
#' @param delete_tar logical whether to delete the intermediate 
#' compressed datasets.
#' See \code{\link{preprocessed_10x_data}} for details.
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
#' function and \code{\link{preprocessed_10x_data}} function. This function will
#' download the preprocessed 10x datasets from a remote host and load
#' them into R as a SingleCellExperiment object or a list of it if
#' querying multiple datasets. To list all available datasets,
#' simply run \code{preprocessed_10x_data()}
#' This function takes the combination
#' of the parameters of the \code{\link[fishpond]{loadFry}}
#' function and \code{\link{preprocessed_10x_data}} function. For
#' each of the \code{\link[fishpond]{loadFry}} parameters,
#' a list of that parameter can be specified, each for a
#' queried dataset.
#' 
#' @export
#'
#' @return If a single dataset is queried, a SingleCellExperiment
#' object will be returned. If multiple datasets are queried,
#' a list of SingleCellExperiment objects will be returned, each
#' for a queried dataset.
#'
#' @examples
#' 
#' \dontrun{
#' library(roe)
#' # run the function
#' # The four way to define output_format
#' # are the same.
#' load10x(dataset_ids = c(1, 2),
#'         output_dir = "10x_datasets",
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

load10x <- function(dataset_ids,
                    output_dir = "preprocessed_10x_data",
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
    # output_format <- list(list(counts = c("U", "S")),
    #                         list(counts = c("U", "S")),
    #                         list(counts = c("U", "S")))
    if (is.list(output_format[[1]])) {
        # user defines outputFormat for each dataset
        # check whether length matches
        if (length(output_format) != nd) {
            stop("The providing output_format list has different length with dataset_ids, cannot proceed")
        }
    } else if (is.list(output_format)) {
        # if user provides a list, it can
        # either be a customized outputFormat
        # or a list of pre-defined outputFormats
        # if customized, then repeat it for each
        # queried dataset
        # output_format <- list(counts = c("U", "S"))
        if (sum(output_format[[1]] %in% c("U", "S", "A")) > 0) {
            # if user provides a customized outputFormat,
            # repeat it for each queried dataset
            # output_format <- rep(output_format, nd)
        } else {
            # if user set a outputFormat for eahc queried dataset,
            # then check the length
            # output_format <- list("scRNA", "scRNA", "scRNA")
            if (length(output_format) != nd) {
                stop("The providing output_format list has different length with dataset_ids, cannot proceed")
            }
        }
    } else {
        # output_format <- "scRNA"
        # The last valid situation is that the user provides a sinlg pre-defined format
        # then we repeat it for every dataset
        output_format <- rep(output_format, nd)
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
        nonzero <- rep(nonzero, nd)
    }

    .say(quiet, "Downloading datasets")

    # download the datsets
    dataset_paths <- preprocessed_10x_data(dataset_ids = dataset_ids,
                            output_dir = output_dir,
                            force = force,
                            delete_tar =  delete_tar,
                            quiet = quiet
                            )

    sce_list <- list()
    # process them using user output
    for (dataset_id in seq(nd)) {
        .say(quiet, "Loading dataset ", dataset_ids[dataset_id])

        dataset_path_ds <- dataset_paths[dataset_id]
        output_format_ds <- output_format[[dataset_id]]
        nonzero_ds <- nonzero[[dataset_id]]
        sce_list[[dataset_id]] <-
                            fishpond::loadFry(fryDir = dataset_path_ds,
                                                outputFormat = output_format_ds,
                                                nonzero = nonzero_ds,
                                                quiet = quiet
                        )
    }
}

