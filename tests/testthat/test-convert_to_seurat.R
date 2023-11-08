data(mbrain_raw)

spatial_dir <- system.file(file.path("extdata",
                                     "V1_Adult_Mouse_Brain_spatial"),
                           package = "SpotClean")
mbrain_slide_info <- read10xSlide(tissue_csv_file=file.path(spatial_dir,
                                       "tissue_positions_list.csv"),
             tissue_img_file = file.path(spatial_dir,
                                       "tissue_lowres_image.png"),
             scale_factor_file = file.path(spatial_dir,
                                       "scalefactors_json.json"))

mbrain_obj <- createSlide(mbrain_raw, mbrain_slide_info)

test_that("Invalid image directory",{
    expect_error(convertToSeurat(mbrain_obj, "foo", "raw"),
                 "unable to open foo/tissue_lowres_image.png")
})

seurat_obj <- convertToSeurat(mbrain_obj,image_dir = spatial_dir,
                              filter_matrix = FALSE)

test_that("Consistent count matrix",{
    if(as.integer(gsub("\\<(\\d+)\\.\\d+\\.\\d+", "\\1", 
                       seurat_obj@version))>=5){
        seurat_counts <- seurat_obj@assays$Spatial$counts
    }else{
        seurat_counts <- seurat_obj@assays$Spatial@counts
    }
    
    expect_identical(mbrain_obj@assays@data$raw,
                     seurat_counts)
})

test_that("Consistent slide metadata",{
    expect_identical(factor(seurat_obj@images$slice1@coordinates$tissue),
    mbrain_obj@metadata$slide$tissue)
    expect_identical(rownames(seurat_obj@images$slice1@coordinates),
                     mbrain_obj@metadata$slide$barcode)
    expect_identical(seurat_obj@images$slice1@coordinates[,c("row","col")],
                     mbrain_obj@metadata$slide[,c("row","col")])
})

seurat_obj_f <- convertToSeurat(mbrain_obj,image_dir = spatial_dir,
                              filter_matrix = TRUE)

test_that("Filter background spots",{
    expect_equal(dim(seurat_obj_f), c(100, 2702))
    expect_identical(
        sort(colnames(seurat_obj_f)),
        sort(mbrain_slide_info$slide$barcode[mbrain_slide_info$slide$tissue==1])
        )
})
