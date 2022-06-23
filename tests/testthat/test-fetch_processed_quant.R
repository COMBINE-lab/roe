suppressPackageStartupMessages({
    library(roe)
    library(SingleCellExperiment)
})
test_that("fetch processed quant function works", {
    output_dir <- tempdir()
    ds <- get_available_dataset_df()[22, ]
    pq <- fetch_processed_quant(22, fetch_dir = output_dir)[[1]]
    tar_dir <- file.path(output_dir, "quant_tar")
    expected_tar_path <- file.path(tar_dir, paste0(pq@dataset_id, ".tar"))
    quant_dir <- output_dir
    expected_quant_path <- file.path(quant_dir,
        pq@dataset_id,
        paste0(pq@fastq_MD5sum,
            "_fry_unfilt_quant_usa_cr-like"
        )
    )
    expect_equal(pq@dataset_id,  ds$dataset_id)
    expect_equal(pq@chemistry,  ds$chemistry)
    expect_equal(pq@reference,  ds$reference)
    expect_equal(pq@dataset_name,  ds$dataset_name)
    expect_equal(pq@dataset_url,  ds$dataset_url)
    expect_equal(pq@fastq_url,  ds$fastq_url)
    expect_equal(pq@fastq_MD5sum,  ds$fastq_MD5sum)
    expect_equal(pq@feature_barcode_csv_url, character(0))
    expect_equal(pq@multiplexing_library_csv_url, character(0))
    expect_equal(pq@quant_tar_url,  ds$quant_tar_url)
    expect_equal(file.exists(expected_tar_path), TRUE)
    expect_equal(file.exists(expected_quant_path), TRUE)
    expect_s4_class(pq@sce, "SingleCellExperiment")

})