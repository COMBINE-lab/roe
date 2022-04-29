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
#' @param keep_tar logical whether to delete the intermediate 
#' compressed datasets.
#' See \code{\link{fetch_processed_quant}} for details.
#' @param output_format can be \emph{either} a valid
#' \code{\link[fishpond]{loadFry}} \code{outputFormat}
#' or a list of it. The \code{names} of the list must
#' match the \code{dataset_ids}.
#' @param nonzero Similar with \code{outputFormat}, it can be
#' either a valid \code{\link[fishpond]{loadFry}} \code{nonzero}
#' parameter or a list of it. If a list, The \code{names} of the 
#' list must match the \code{dataset_ids}.
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
#' a data frame containing the information of available datasets
#' will be returned; otherwise, a list of lists, where each list
#' stores the information of one fetched dataset. The `quant_dir`
#' field represent the path to the quantification result of the 
#' fetched dataset; The `sce` field represent the SingleCellExperiment
#' of this dataset.
#'
#' @examples
#' 
#' \dontrun{
#' library(roe)
#' # run the function
#' # The four ways to define output_format
#' # are the same.
#' load_processed_quant(dataset_ids = c(1, 2),
#'         fetch_dir = "processed_quant",
#'         force = FALSE,
#'         keep_tar = TRUE,
#'         output_format = "scRNA",
#' #         output_format = list("scRNA", "scRNA"),
#' #         output_format = list("1" = list(counts = c("S", "A")),
#' #                               "2" = list(counts = c("S", "A"))
#' #                              ),
#' #         output_format = list("counts" = c("S", "A")),
#'         nonzero = FALSE,
#'         quiet = FALSE
#' )
#' }
#' 

load_processed_quant <- function(dataset_ids = c(),
                    fetch_dir = "processed_quant",
                    force = FALSE,
                    keep_tar = TRUE,
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

    dataset_ids = check_dataset_ids(dataset_ids)
    # check whether there is any dataset id left
    if (length(dataset_ids) == 0) {
        stop("No valid dataset id found, can not proceed")
    }
    
    # check output_format
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
    
    pq_list <- list()
    # folder for (temporarily) storing tar files.
    tar_dir <- file.path(fetch_dir, "datasets_tar")
    dir.create(tar_dir, recursive = TRUE,
               showWarnings = FALSE)
    
    # process them using user output
    for (dataset_id in dataset_ids) {
        
        output_format_ds <- output_format[[dataset_id]]
        nonzero_ds <- nonzero[[dataset_id]]

        processed_dataset = FDL(dataset_id,
                                tar_dir=tar_dir,
                                quant_dir=fetch_dir,
                                output_format=output_format_ds,
                                nonzero=nonzero_ds,
                                force=force, 
                                quiet=quiet)

        
        
        # reset tar_path if needed
        if (!keep_tar) {
            processed_quant$tar_path = NULL
        }
        
        # append to list
        pq_list[[as.character(dataset_id)]] = processed_quant
    }
    
    if (!keep_tar) {
        .say(quiet,
             "Removing downloaded tar files in directory:\n",
             "  ", tar_dir)
        unlink(tar_dir,  recursive = TRUE, force = TRUE)
    }
    
    .say(quiet, "Done")
    return(pq_list)
}
