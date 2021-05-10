# Simulate 10x output files

data(mbrain_raw)
data(mbrain_slide_info)

data_dir <- file.path(tempdir(),"sim_example")
if(dir.exists(data_dir)){
    file.remove(list.files(data_dir,full.names = TRUE))
}else{
    dir.create(data_dir)
}

matrix_dir <- file.path(data_dir,"matrix.mtx")
barcode_dir <- gzfile(file.path(data_dir, "barcodes.tsv.gz"), open="wb")
gene_dir <- gzfile(file.path(data_dir, "features.tsv.gz"), open="wb")

# For simplicity, use gene names to generate gene IDs to fit the format.
gene_name <- rownames(mbrain_raw)
gene_id <- paste0("ENSG_fake_",gene_name)
barcode_id <- colnames(mbrain_raw)

Matrix::writeMM(mbrain_raw,file = matrix_dir)
write(barcode_id,file = barcode_dir)
write.table(cbind(gene_id,gene_name,"type"),file = gene_dir,
    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
R.utils::gzip(matrix_dir)
close(barcode_dir)
close(gene_dir)

# test reading function

test_that("Directory/files existence", {
    expect_error(Read10xRaw(paste0(data_dir,"/foo")),
                 "Directory does not exist")
    expect_error(Read10xRaw(tempdir()), "No 10x output file detected")
})

test_that("Data loading", {
    expect_identical(Read10xRaw(data_dir), mbrain_raw)
})


test_that("Metadata loading", {
    expect_true(is.list(Read10xRaw(data_dir, meta = TRUE)))
})
