.say <- function(quiet, ...) {
    if (!quiet) {
        message(...)
    }
}

check_dataset_ids <- function(dataset_ids) {
  # now check the validity of dataset_ids
  if (is.numeric(dataset_ids)) {
    valid_ids <- (dataset_ids >= 1) &
      (dataset_ids <= nrow(available_datasets))
    if (any(!valid_ids)) {
      message("Found invalid dataset id: '",
              paste0(dataset_ids[!valid_ids], collapse = "' '"),
              "', ignored")
      dataset_ids <- dataset_ids[valid_ids]
    }
  } else {
    stop("dataset ids must be integer, can not proceed")
  }
  return(dataset_ids)
}