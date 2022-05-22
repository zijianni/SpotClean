#' @title Convert slide object to Seurat object
#'
#' @description This function converts our slide object of class
#' \code{SummarizedExperiment} to Seurat object of class \code{Seurat} so that
#' users can directly proceed with Seurat spatial analyses pipelines.
#' Built based on Seurat's \code{Load10X_Spatial()}.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{CreateSlide()}.
#'
#' @param image_dir (chr) Path to directory with 10X Genomics visium image data;
#' should include files \code{tissue_lowres_iamge.png},
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
#' \dontrun{
#' data(mbrain_raw)
#' data(mbrain_slide_info)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' example_image_dir <- "path/to/image/dir"
#' seurat_obj <- convertToSeurat(mbrain_obj, example_image_dir, "raw")
#' str(seurat_obj)
#' }
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
