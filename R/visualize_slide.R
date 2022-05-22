#' @title Visualize the Visium slide image
#'
#' @description Generate and visualize the tissue image as a ggplot2 object.
#' Users can manually add and modify layers (e.g. title, axis)
#' following ggplot2's syntax.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{createSlide()}.
#'
#' @param title (chr) Title of the plot. Default: \code{""}
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#'
#' data(mbrain_raw)
#' data(mbrain_slide_info)
#' mbrain_obj <- createSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' gp <- visualizeSlide(mbrain_obj)
#' plot(gp)


#' @importFrom ggplot2 ggplot aes coord_cartesian xlim ylim xlab ylab
#' @importFrom ggplot2 ggtitle guides guide_legend theme_set theme_bw
#' @importFrom ggplot2 theme element_blank element_line
#' @importFrom tibble tibble
#' @importFrom S4Vectors metadata
#'
#' @export

visualizeSlide <- function(slide_obj, title=""){

    # junk code... get rid of R CMD check notes
    imagerow <- imagecol <- NULL

    slide <- metadata(slide_obj)$slide
    grob <- metadata(slide_obj)$grob

    if(is.null(grob)){
        stop("No valid image information. ",
        "Check path to image file in read10xSlide().")
    }

    gp <- ggplot(slide, aes(x = imagecol, y = imagerow)) +
        .geom_spatial(data = tibble(grob=list(grob)),
                      aes(grob = grob),
                      x = 0.5,
                      y = 0.5) +
        coord_cartesian(expand = FALSE) +
        xlim(0, max(slide$width)) +
        ylim(max(slide$height), 0) +
        xlab("") +
        ylab("") +
        ggtitle(title)+
        guides(fill = guide_legend(override.aes = list(size = 3))) +
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


#' @import grid
#' @importFrom ggplot2 ggproto Geom ggproto_parent layer
# This function is developed based on 10x's secondary analysis pipeline
# https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/rkit.

.geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {

    GeomCustom <- ggproto(
        "GeomCustom",
        Geom,
        setup_data = function(self, data, params) {
            data <- ggproto_parent(Geom, self)$setup_data(data, params)
            data
        },

        draw_group = function(data, panel_scales, coord) {
            vp <- viewport(x=data$x, y=data$y)
            g <- editGrob(data$grob[[1]], vp=vp)
            .ggname(".geom_spatial", g)
        },

        required_aes = c("grob","x","y")

    )

    layer(
        geom = GeomCustom,
        mapping = mapping,
        data = data,
        stat = stat,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}


# Name ggplot grid object, from ggplot2/R/utilities-grid.r
# Convenience function to name grid objects

.ggname <- function(prefix, grob) {
    grob$name <- grobName(grob, prefix)
    grob
}
