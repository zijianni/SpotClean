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
#' @param image_dir Path to directory with 10X Genomics visium image data;
#' should include files \code{tissue_lowres_iamge.png},
#' \code{scalefactors_json.json} and \code{tissue_positions_list.csv}.
#'
#'
#' @param slice Name for the stored image of the tissue slice.
#' Default: "slice1".
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
#' example_image_dir <- 'path/to/image/dir'
#' seurat_obj <- ConvertToSeurat(mbrain_obj, example_image_dir, "raw")
#' str(seurat_obj)
#' }
#'
#' @importFrom Seurat Read10X_Image Cells CreateSeuratObject DefaultAssay<-
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @export

ConvertToSeurat <- function(slide_obj, image_dir,
                            slice="slice1"){


    object <- CreateSeuratObject(assay(slide_obj),
                                 assay = "Spatial")

    filter_matrix <- ncol(object)!=nrow(metadata(slide_obj)$slide)

    image <- Read10X_Image(image.dir = image_dir,
                           filter.matrix = filter_matrix)
    image <- image[Cells(x = object)]
    DefaultAssay(object = image) <- "Spatial"
    object[[slice]] <- image
    return(object)
}
