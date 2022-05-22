#' @title Visualize spot labels on the 2D slide
#'
#' @description Generate and visualize spot labels on the 2D slide
#' as a ggplot2 object. Spot labels can be given via parameter \code{label} as
#' external input or one column in the \code{slide} slot of the slide
#' object's metadata via \code{label_col}. Exactly one of \code{label}
#' and \code{label_col} must be specified.
#'
#' @param object A slide object created or
#' inherited from \code{createSlide()}, or a \code{data.frame} of slide
#' information with columns: barcodes, tissue, imagerow, imagecol, etc.
#'
#' @param ... Arguments passed to other methods
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#'
#' data(mbrain_raw)
#' data(mbrain_slide_info)
#' mbrain_obj <- createSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' gp <- visualizeLabel(mbrain_obj, label="tissue",
#'                      title="mbrain", legend_title="tissue or background")
#' plot(gp)
#'
#' @rdname visualizeLabel
#'
#' @export

visualizeLabel <- function(object, ...) {
    UseMethod(generic = "visualizeLabel", object = object)
}


#' @param label (categorical vector or chr) Either a vector of labels for all
#' spots, or the column name in the \code{slide} slot of the slide
#' object's metadata. In the former case, the order of values in the vector
#' should match the spot barcodes in the slide object. In the latter case,
#' The column should be categorical with reasonably number of unique labels.
#' Default: \code{"tissue"}
#'
#' @param subset_barcodes (vector of chr) A subset of spot barcodes to plot.
#' By default it plots all spots in the slide object. This can be useful
#' when only plotting tissue spots or specific tissue types or regions.
#' Default: \code{NULL}
#'
#' @param title (chr) Title of the plot. Default: \code{""}
#'
#' @param legend_title (chr) Title of the legend. Under default,
#' use \code{label} as legend title. Default: \code{NULL}
#'


#' @importFrom methods as
#'
#' @method visualizeLabel default
#' @rdname visualizeLabel
#'
#' @export

visualizeLabel.default <- function(object, label="tissue",
                                   subset_barcodes=NULL, title="",
                                   legend_title=NULL, ...){

    if (!inherits(x = object, "data.frame")) {
        object <- as(object = object, Class = "data.frame")
    }

    slide <- object

    gp <- .visualize_label(slide, label,
                           subset_barcodes, title,
                           legend_title, ...)
    return(gp)

}

#' @method visualizeLabel SummarizedExperiment
#' @rdname visualizeLabel
#' @importFrom S4Vectors metadata
#' @export
#'
visualizeLabel.SummarizedExperiment <- function(object, label="tissue",
                                                subset_barcodes=NULL, title="",
                                                legend_title=NULL, ...){

    slide <- metadata(object)$slide

    gp <- .visualize_label(slide, label,
                           subset_barcodes, title,
                           legend_title, ...)
    return(gp)

}

#' @importFrom ggplot2 ggplot aes geom_point coord_cartesian xlim ylim
#' @importFrom ggplot2 xlab ylab ggtitle labs guides guide_legend
#' @importFrom ggplot2 theme_set theme_bw theme element_blank element_line
#' @importFrom dplyr filter
#'
.visualize_label <- function(slide, label="tissue",
                             subset_barcodes=NULL, title="",
                             legend_title=NULL, ...){

    # junk code... get rid of R CMD check notes
    imagerow <- imagecol <- barcode <- NULL

    # manipulate label to plot
    if(length(label)==1){

        if(!label%in%colnames(slide)){
            stop("Label name does not exist in slide metadata.")
        }
        if(is.null(legend_title)){
            legend_title <- label
        }
        label <- slide[,label]

    }else if(length(label)==nrow(slide)){

        label <- as.character(label)

    }else{
        stop("Invalid label input.")
    }

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
