library(roe)
library(SingleCellExperiment)
test_that("ProcesedQuant class initialization works", {
    ds <- get_available_dataset_df()[22, ]
    pq <- ProcessedQuant(22)
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
    expect_equal(pq@quant_path, character(0))
    expect_equal(pq@tar_path, character(0))
    expect_equal(pq@sce, SingleCellExperiment())
    }
)

test_that("ProcessedQuant class can fetch, decompress and load a dataset", {
    output_dir <- tempdir()
    pq <- ProcessedQuant(22)

    # test error
    expect_error(decompress_quant(pq, quiet = FALSE))
    expect_error(load_quant(pq, quiet = FALSE))


    # test fetch_quant
    tar_dir <- file.path(output_dir, "quant_tar")
    tar_name <- "test_ds22_quant_tar"
    pq <- fetch_quant(pq, tar_dir = tar_dir, tar_name = tar_name, quiet = FALSE)
    expected_tar_path <- file.path(tar_dir, paste0(tar_name, ".tar"))
    expect_equal(file.exists(expected_tar_path), TRUE)

    ## test force
    pq <- fetch_quant(pq,
        tar_dir = tar_dir,
        tar_name = tar_name,
        force = TRUE,
        quiet = FALSE
    )
    expect_equal(file.exists(expected_tar_path), TRUE)

    # test error
    expect_error(load_quant(pq))

    # test decompress_quant
    quant_path_name <- "test_ds22_quant_dir"
    quant_dir <- output_dir
    pq <- decompress_quant(pq,
        quant_dir = quant_dir,
        quant_path_name = "test_ds22_quant_dir",
        quiet = FALSE
    )

    expected_quant_path <- file.path(quant_dir,
        quant_path_name,
        paste0(pq@fastq_MD5sum,
            "_fry_unfilt_quant_usa_cr-like"
        )
    )
    expect_equal(file.exists(expected_quant_path), TRUE)

    ## test force
    pq <- decompress_quant(pq,
        quant_dir = quant_dir,
        quant_path_name = "test_ds22_quant_dir",
        force = TRUE,
        quiet = FALSE
    )
    expect_equal(file.exists(expected_quant_path), TRUE)

    # test load_quant
    pq <- load_quant(pq, quiet = FALSE)
    expect_s4_class(pq@sce, "SingleCellExperiment")

    # test a few more output_formats
    output_format <- "raw"
    pq <- load_quant(pq,
        output_format = output_format,
        quiet = FALSE
    )
    expected_assay_name <- c("spliced", "unspliced", "ambiguous")
    expect_equal(sort(expected_assay_name),
        sort(names(pq@sce@assays))
    )

    output_format <- list(testA = c("S", "A"),
        testB = c("U")
    )
    pq <- load_quant(pq,
        output_format = output_format,
        quiet = FALSE)
    expected_assay_name <- c("testA", "testB")
    expect_equal(sort(expected_assay_name),
        sort(names(pq@sce@assays))
    )


})
