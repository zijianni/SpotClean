#' @title Visualize spot labels on the 2D slide
#'
#' @description Generate and visualize spot labels on the 2D slide
#' as a ggplot2 object. Spot labels can be given via parameter \code{label} as
#' external input or one column in the \code{slide} slot of the slide
#' object's metadata via \code{label_col}. Exactly one of \code{label}
#' and \code{label_col} must be specified.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{CreateSlide()}.
#'
#' @param label (categorical vector or chr) Either a vector of labels for all
#' spots, or the column name in the \code{slide} slot of the slide
#' object's metadata. In the former case, the order of values in the vector
#' should match the spot barcodes in the slide object. In the latter case,
#' The column should be categorical with reasonably number of unique labels.
#' Default: \code{"tissue"}.
#'
#' @param subset_barcodes (vector of chr) A subset of spot barcodes to plot.
#' By default it plots all spots in the slide object. This can be useful
#' when only plotting tissue spots or specific tissue types or regions.
#' Default: \code{NULL}.
#'
#' @param title (chr) Title of the plot. Default: \code{""}.
#'
#' @param legend_title (chr) Title of the legend. Default: \code{"Label"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#'
#' data(MbrainSmall)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' gp <- VisualizeLabel(mbrain_obj, label="tissue",
#'                      title="mbrain", legend_title="tissue or background")
#' plot(gp)


#' @import ggplot2
#' @import dplyr
#' @importFrom SummarizedExperiment metadata
#'
#' @export

VisualizeLabel <- function(slide_obj, label="tissue",
                           subset_barcodes=NULL, title="",
                           legend_title="Label"){

    # manipulate label to plot
    if(length(label)==1){

        if(!label%in%colnames(metadata(slide_obj)$slide)){
            stop("Label does not exist in slide metadata.")
        }

        label <- metadata(slide_obj)$slide[,label]

    }else if(length(label)==ncol(slide_obj)){

        label <- as.character(label)

    }else{
        stop("Invalid label input.")
    }

    slide <- metadata(slide_obj)$slide
    slide$label <- label

    # subsetting barcodes
    if(!is.null(subset_barcodes)){
        slide <- filter(slide, barcode%in%subset_barcodes)
    }

    # plot
    gp <- ggplot(slide, aes(x=imagecol,y=imagerow, fill=factor(label))) +
        geom_point(shape = 21, size = 1.75, color="white")+
        coord_cartesian(expand=FALSE)+
        xlim(0,max(slide$width))+
        ylim(max(slide$height),0)+
        xlab("") +
        ylab("") +
        ggtitle(title)+
        labs(fill = legend_title)+
        guides(fill = guide_legend(override.aes = list(size=3)))+
        theme_set(theme_bw(base_size = 10))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text = element_blank(),
              axis.ticks = element_blank())

    return(gp)

}
