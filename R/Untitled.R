methods::setClass("ProcessedQuant", methods::representation(dataset_id = "numeric",
                                          chemistry = "character",
                                          reference = "character",
                                          dataset_name = "character",
                                          link = "character",
                                          data_url = "character",
                                          MD5	= "character",
                                          feature_barcode	= "character",
                                          library_csv	= "character",
                                          quant_link = "character",
                                          quant_path = "character",
                                          tar_path = "character",
                                          sce = "SingleCellExperiment"))


load("R/sysdata.rda")
dataset_id = 2
dataset_info <- available_datasets[dataset_id, ]
library(SingleCellExperiment)
pq = methods::new("ProcessedQuant",
    dataset_id = dataset_info$"dataset_id",
    chemistry = dataset_info$"chemistry",
    reference = dataset_info$"reference",
    dataset_name = dataset_info$"dataset_name",
    link = dataset_info$"link",
    data_url = dataset_info$"data_url",
    MD5	= dataset_info$"MD5",
    sce = SingleCellExperiment::SingleCellExperiment())

if (!is.na(dataset_info$feature_barcode)) {
    pq@feature_barcode = dataset_info$feature_barcode
}
if (!is.na(dataset_info$library_csv)) {
    pq@library_csv = dataset_info$library_csv
}

setMethod("sce", "ProcessedQuant", function(x) x@sce)





# new("ProcessedQuant",
#     dataset_id = dataset_info$"dataset_id",
#     chemistry = dataset_info$"chemistry",
#     reference = dataset_info$"reference",
#     dataset_name = dataset_info$"dataset_name",
#     link = dataset_info$"link",
#     data_url = dataset_info$"data_url",
#     MD5	= dataset_info$"MD5",
#     feature_barcode	= dataset_info$"feature_barcode",
#     library_csv	= dataset_info$"library_csv",
#     quant_link = dataset_info$"quant_link",
#     quant_path = NULL,
#     tar_path = NULL,
#     sce = NULL)
