#' @title Visualize spot values on the 2D slide
#'
#' @description Generate and visualize spot labels on the 2D slide
#' as a ggplot2 object. Spot labels can be given via parameter \code{label} as
#' external input or one column in the \code{slide} slot of the slide
#' object's metadata via \code{label_col}. Exactly one of \code{label}
#' and \code{label_col} must be specified.
#'
#' @param object A slide object created or
#' inherited from \code{CreateSlide()}, or a \code{data.frame} of slide
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
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' gp <- VisualizeHeatmap(mbrain_obj, "Bc1",
#'                      title="mbrain", legend_title="Bc1 expression")
#' plot(gp)


#' @rdname VisualizeHeatmap
#'
#' @export

VisualizeHeatmap <- function(object, ...) {
    UseMethod(generic = "VisualizeHeatmap", object = object)
}


#' @param value (numeric vector or chr) Either a vector of numeric values for
#' all spots, or the name of one gene in \code{exp_matrix} or in the expression
#' matrix in the slide object.
#' In the former case, the order of values in the
#' vector should match the spot barcodes in the slide object.
#'
#' @param exp_matrix (matrix of num) When \code{value} is a single character
#' string, will search the matching gene in \code{exp_matrix} and plot the gene
#' expression. Default: \code{NULL}.
#'
#' @param subset_barcodes (vector of chr) A subset of spot barcodes to plot.
#' By default it plots all spots in the slide object. This can be useful
#' when only plotting tissue spots or specific tissue types or regions.
#' Default: \code{NULL}.
#'
#' @param logged (logical) Specify if the color scale is log1p transformed.
#' Default: \code{TRUE}.
#'
#' @param viridis (logical) If true, color scale uses viridis.
#' Otherwise, use rainbow. Default: \code{TRUE}.
#'
#' @param legend_range (length 2 vector of num) Custom legend range of the value.
#' By default uses the range of the plotted values.
#' Default: \code{NULL}.
#'
#' @param title (chr) Title of the plot. Default: \code{""}.
#'
#' @param legend_title (chr) Title of the legend. Under default,
#' use \code{value} as legend title. Default: \code{NULL}.
#'
#' @importFrom methods as
#'
#' @method VisualizeHeatmap default
#' @rdname VisualizeHeatmap
#'
#' @export

VisualizeHeatmap.default <- function(object, value, exp_matrix=NULL,
                             subset_barcodes=NULL,
                             logged=TRUE, viridis=TRUE,
                             legend_range=NULL,
                             title="", legend_title=NULL,
                             ...){

    if (!inherits(x = object, "data.frame")) {
        object <- as(object = object, Class = "data.frame")
    }

    slide <- object

    gp <- .visualize_heatmap(slide, value, exp_matrix,
                             subset_barcodes,
                             logged, viridis,
                             legend_range,
                             title, legend_title,
                             ...)
    return(gp)
}


#' @method VisualizeHeatmap SummarizedExperiment
#' @rdname VisualizeHeatmap
#' @importFrom SummarizedExperiment assay
#' @importMethodsFrom S4Vectors metadata
#' @export
#'
VisualizeHeatmap.SummarizedExperiment <- function(object, value,
                                     subset_barcodes=NULL,
                                     logged=TRUE, viridis=TRUE,
                                     legend_range=NULL,
                                     title="", legend_title=NULL,
                                     ...){

    slide <- metadata(object)$slide
    exp_matrix <- assay(object)

    gp <- .visualize_heatmap(slide, value, exp_matrix,
                       subset_barcodes,
                       logged, viridis,
                       legend_range,
                       title, legend_title,
                       ...)
    return(gp)
}


#' @import ggplot2
#' @importFrom dplyr filter
#'
.visualize_heatmap <- function(slide, value, exp_matrix=NULL,
                               subset_barcodes=NULL,
                               logged=TRUE, viridis=TRUE,
                               legend_range=NULL,
                               title="", legend_title=NULL,
                               ...){
    # junk code... get rid of R CMD check notes
    imagerow <- imagecol <- barcode <- NULL


    # manipulate value to plot
    if(length(value)==1 & is.character(value)){
        if(is.null(legend_title)){
            legend_title <- value
        }

        # value to plot is in slide dataframe
        if(value%in%colnames(slide)){
            slide$value <- slide[,value]
        }else{
            # value to plot is in given matrix
            if(is.null(exp_matrix)){
                stop("You must provide an input expression matrix to plot ",
                     value," expressions.")
            }
            if(!value%in%rownames(exp_matrix)){
                stop("Specified gene does not exist in the expression matrix.")
            }

            # if expression matrix does not match slide info
            shared_bcs <- intersect(colnames(exp_matrix), slide$barcode)
            if(length(shared_bcs)==0){
                stop("Barcodes in input matrix do not ",
                     "match any barcodes in slide.")
            }
            missed_bcs <- setdiff(slide$barcode,shared_bcs)

            value <- c(exp_matrix[value,shared_bcs],rep(NA, length(missed_bcs)))
            names(value) <- c(shared_bcs, missed_bcs)
            slide$value <- value[slide$barcode]
        }

    }else if(length(value)==nrow(slide)){
        # values to plot is directly given
        slide$value <- as.numeric(value)
    }else{
        stop("Invalid value input.")
    }

    # subsetting barcodes
    if(!is.null(subset_barcodes)){
        slide <- filter(slide, barcode%in%subset_barcodes)
    }

    # setup legend breaks
    if(is.null(legend_range)){
        legend_range <- c(min(0,min(slide$value, na.rm = TRUE)),
                          max(slide$value, na.rm = TRUE))
    }
    if(logged){
        legend_breaks <- floor(expm1( log1p(min(legend_range))+
                                          diff(log1p(legend_range))/4*0:4 ))
    }else{
        legend_breaks <- floor(min(legend_range)+diff(legend_range)/4*0:4)
    }

    if(min(legend_range)==0){
        slide_show <- filter(slide,value>0)
        slide_hide <- filter(slide, value==0)
    }else{
        slide_show <- filter(slide,
                             value>=min(legend_range),value<=max(legend_range))
        slide_hide <- filter(slide,
                             value<min(legend_range) | value>max(legend_range))
        if(nrow(slide_hide)>0){
            warning(nrow(slide_hide), " spots outside specified range.\n")
        }

    }

    # plot
    gp <- ggplot(slide_show, aes(x = imagecol, y = imagerow, fill = value)) +
        geom_point(
            shape = 21,
            colour = "white",
            size = 1.75) +
        geom_point(data=slide_hide,
                   shape = 21,
                   colour = "white",
                   fill="#d6d6d6",
                   size = 1.75,
                   alpha=0.3
        ) +
        coord_cartesian(expand = FALSE) +
        .scale_fill_fun(viridis=viridis,
                        trans = ifelse(logged,"log1p","identity"),
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

# Color scales: viridis v.s. rainbow
#' @importFrom viridis scale_fill_viridis
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
.scale_fill_fun <- function(viridis=TRUE, ...){
    if(viridis){
        return(scale_fill_viridis(...))
    }else{
        myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
        return(
            scale_fill_gradientn(colours = myPalette(100), ...)
        )
    }
}
