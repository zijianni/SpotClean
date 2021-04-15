#' @title Visualize spot values on the 2D slide
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
#' @param value (numeric vector or chr) Either a vector of numeric values for all
#' spots, or the name of one gene in an assay specified in \code{assay}.
#' In the former case, the order of values in the vector
#' should match the spot barcodes in the slide object.
#'
#' @param subset_barcodes (vector of chr) A subset of spot barcodes to plot.
#' By default it plots all spots in the slide object. This can be useful
#' when only plotting tissue spots or specific tissue types or regions.
#' Default: \code{NULL}.
#'
#' @param assay_name (chr) When plotting gene expressions, specify the assay
#' to use in the slide object. Available assay names can be found using \
#' \code{assayNames(<slide_obj>)}.
#' Default: \code{"raw"}.
#'
#' @param logged (logical) Specify if the color scale is log1p transformed.
#' Default: \code{TRUE}.
#'
#' @param legend_range (length 2 vector of num) Custom legend range of the value.
#' By default uses the range of the plotted values.
#' Default: \code{NULL}.
#'
#' @param title (chr) Title of the plot. Default: \code{""}.
#'
#' @param legend_title (chr) Title of the legend. Default: \code{"Value"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#'
#' data(MbrainSmall)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' gp <- VisualizeHeatmap(mbrain_obj, "Bc1",
#'                      title="mbrain", legend_title="Bc1 expression")
#' plot(gp)


#' @import ggplot2
#' @import viridis
#' @import dplyr
#' @importFrom SummarizedExperiment assays metadata
#'
#' @export

VisualizeHeatmap <- function(slide_obj, value,
                             subset_barcodes=NULL, assay_name="raw",
                             logged=TRUE, legend_range=NULL,
                             title="", legend_title="Value"){

    # manipulate value to plot
    if(length(value)==1){
        count_mat <- assays(slide_obj)[[assay_name]]

        if(!value%in%rownames(count_mat)){
            stop("Specified gene does not exist in the expression matrix.")
        }

        value <- count_mat[value,]

    }else if(length(value)==ncol(slide_obj)){

        value <- as.numeric(value)

    }else{
        stop("Invalid value input.")
    }

    slide <- metadata(slide_obj)$slide
    slide$value <- value

    # subsetting barcodes
    if(!is.null(subset_barcodes)){
        slide <- filter(slide, barcode%in%subset_barcodes)
    }

    # setup legend breaks
    if(is.null(legend_range)){
        legend_range <- range(slide$value)
    }
    if(logged){
        legend_breaks <- floor(expm1( log1p(min(legend_range))+
                                          diff(log1p(legend_range))/4*0:4 ))
    }else{
        legend_breaks <- floor(min(legend_range))+diff(legend_range)/4*0:4
    }

    # plot
    gp <- ggplot(slide %>% filter(value>0), aes(x = imagecol, y = imagerow, fill = value)) +
        geom_point(
            shape = 21,
            colour = "white",
            size = 1.75) +
        geom_point(data=slide %>% filter(value==0),
               shape = 21,
               colour = "white",
               fill="#d6d6d6",
               size = 1.75,
               alpha=0.3
        ) +
        coord_cartesian(expand = FALSE) +
        scale_fill_viridis(trans = ifelse(logged,"log1p","identity"),
                           breaks=legend_breaks, limits=legend_range) +
        xlim(0, max(slide$width)) +
        ylim(max(slide$height), 0) +
        xlab("") +
        ylab("") +
        ggtitle(title) +
        labs(fill = legend_title) +
        theme_set(theme_bw(base_size = 10)) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank()
        )

    return(gp)
}
