.say <- function(quiet, ...) {
    if (!quiet) {
        message(...)
    }
}

.get_dataset_info_list<- function(available_datasets, dataset_id) {    
  processed_quant = list(dataset_id = available_datasets[dataset_id, "dataset_id"],
                         chemistry = available_datasets[dataset_id, "chemistry"],
                         reference = available_datasets[dataset_id, "reference"],
                         dataset_name = available_datasets[dataset_id, "dataset_name"],
                         link = available_datasets[dataset_id, "link"],
                         data_url = available_datasets[dataset_id, "data_url"],
                         MD5	= available_datasets[dataset_id, "MD5"],
                         feature_barcode	= available_datasets[dataset_id, "feature_barcode"],
                         library_csv	= available_datasets[dataset_id, "library_csv"],
                         quant_link = available_datasets[dataset_id, "quant_link"]
  )
  
  processed_quant
}