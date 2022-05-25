#' @title Convert slide object to Seurat object
#'
#' @description This function converts our slide object of class
#' \code{SummarizedExperiment} to Seurat object of class \code{Seurat} so that
#' users can directly proceed with Seurat spatial analyses pipelines.
#' Built based on Seurat's \code{Load10X_Spatial()}.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{createSlide()}.
#'
#' @param image_dir (chr) Path to directory with 10X Genomics visium image data;
#' should include files \code{tissue_lowres_image.png},
#' \code{scalefactors_json.json} and \code{tissue_positions_list.csv}.
#'
#' @param filter_matrix (logical) If \code{TRUE}, only keep spots that have been
#' determined to be over tissue. If \code{slide_obj} only contains tissue spots,
#' \code{filter_matrix} has to be set \code{TRUE}. If \code{slide_obj}
#' contains both tissue and background spots, setting \code{filter_matrix=TRUE}
#' will subset the expression matrix to tissue spots only. Default: \code{TRUE}
#'
#' @param slice (chr) Name for the stored image of the tissue slice.
#' Default: "slice1"
#'
#' @return A Seurat object with spatial information.
#'
#' @examples
#'
#' # load count matrix and slide metadata
#' data(mbrain_raw)
#' spatial_dir <- system.file(file.path("extdata",
#'                                      "V1_Adult_Mouse_Brain_spatial"),
#'                            package = "SpotClean")
#' mbrain_slide_info <- read10xSlide(tissue_csv_file=file.path(spatial_dir,
#'                                        "tissue_positions_list.csv"),
#'              tissue_img_file = file.path(spatial_dir,
#'                                        "tissue_lowres_image.png"),
#'              scale_factor_file = file.path(spatial_dir,
#'                                        "scalefactors_json.json"))
#'
#' # Create slide object
#' mbrain_obj <- createSlide(mbrain_raw,
#'                           mbrain_slide_info)
#'
#' # Convert to Seurat object
#' seurat_obj <- convertToSeurat(mbrain_obj, spatial_dir, "raw")
#' str(seurat_obj)
#'
#' @importFrom Seurat Read10X_Image Cells CreateSeuratObject DefaultAssay<-
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @export

convertToSeurat <- function(slide_obj, image_dir,
                            slice="slice1",
                            filter_matrix=TRUE){

    # create Seurat object
    object <- CreateSeuratObject(assay(slide_obj),
                                 assay = "Spatial")

    # load image and add to Seurat object
    image <- Read10X_Image(image.dir = image_dir,
                           filter.matrix = filter_matrix)
    ts_coord <- GetTissueCoordinates(image)

    if(nrow(ts_coord)>ncol(object)){
        stop("The slide object has fewer spots than the image data.
Consider setting filter_matrix=TRUE?")
    }
    if(nrow(ts_coord)<ncol(object)){
        object <- object[,rownames(ts_coord)]
    }

    image <- image[Cells(x = object)]
    DefaultAssay(object = image) <- "Spatial"
    object[[slice]] <- image

    return(object)
}
