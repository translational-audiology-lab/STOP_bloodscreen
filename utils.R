# ------------------------------------------------------------------------------
# Useful functions
# * irregular_cols
# * geom_split_violin
# * add_n
# ------------------------------------------------------------------------------
# initiated on 2019-12-16
# authors :  Mun-Gwan Hong
# ------------------------------------------------------------------------------

library(tidyverse)

#' Find irregular columns such as too many NAs, just one value
#'
#' @param x a tibble or a data frame
#' @param na_cuff_off a cutoff of the proportion of NAs to find NA enriched
#'   columns
#'
#' @return a list of column names
irregular_cols <- function(x, na_cuff_off = 0.8) {
  out <- list()
  
  # Mostly NA columns
  out$mostly_na <- x %>% 
    select_if(~sum(is.na(.)) > (length(.) * na_cuff_off)) %>% 
    names()
  
  # Just one identical value for all entries of a column
  out$one_value_only <- sapply(names(x), function(ii) {
    col1 <- x[[ii]]
    if(all(!is.na(col1)) & all(col1 == col1[1])) col1[1] else NA
  }) %>% {
    x1 <- .[!sapply(., is.na)]
    tibble(
      Column = names(x1), 
      `Repeated value` = as.character(x1)
    )
  }
  
  # Just one identical value or NA for all entries of a column
  out$one_value_or_na <- names(x)[! names(x) %in% out$one_value_only$Column] %>% 
    sapply(function(ii) {
      col1 <- x[[ii]]
      if(all(is.na(x[[ii]]))) return(NA)
      col1 <- na.omit(col1)
      if(all(col1 == col1[1])) col1[1] else NA
    }) %>% {
      x1 <- .[!sapply(., is.na)]
      tibble(
        Column = names(x1), 
        `Repeated value` = as.character(x1)
      )
    }
  
  # Just one element or less in a group
  out$one_or_less_in_grp <- x %>% 
    select_if(. %>% {is.factor(.) | is.character(.)}) %>% 
    # Remove any factor variable that has a group with only one element
    select_if(function(x1) {
      freqtab <- table(x1)
      # one element in a group
      any(freqtab == 1)
    }) %>% 
    names()
  
  return(out)
}


# -----------------------------------------------------------------------------#
#' Split violin plot
#' 
#' Violin plot showing different distribution on the left and right sides. The
#' code of 'jan-glx' was obtained from StackOverflow (check the reference).
#'
#' @param file a Excel file from the manager
#' 
#' @param mapping,data,stat,position,...,draw_quantiles,trim,scale,na.rm,show.legend,inherit.aes check \code{\link{geom_violin}}
#' 
#' @author jan-glx
#' @references \url{https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2}
#' 
#' @seealso
#' \code{\link{geom_violin}}
#' 
#' @import ggplot2
#' @importFrom scales zero_range
#' @export
# -----------------------------------------------------------------------------#
# created  : 2018-09-12 by Mun-Gwan
# -----------------------------------------------------------------------------#
geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  
  GeomSplitViolin <- ggproto(
    "GeomSplitViolin",
    GeomViolin,
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
      data <-
        transform(
          data,
          xminv = x - violinwidth * (x - xmin),
          xmaxv = x + violinwidth * (xmax - x)
        )
      grp <- data[1, "group"]
      newdata <-
        dplyr::arrange(transform(data, x = if (grp %% 2 == 1)
          xminv
          else
            xmaxv), if (grp %% 2 == 1)
              y
          else-y)
      newdata <-
        rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])
      newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
        round(newdata[1, "x"])
      
      if (length(draw_quantiles) > 0 &
          !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                  1))
        quantiles <-
          ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <-
          data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <-
          rep(1, nrow(quantiles))
        both <-
          cbind(quantiles, aesthetics)
        quantile_grob <-
          GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
                         grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
      }
      else {
        ggplot2:::ggname("geom_split_violin",
                         GeomPolygon$draw_panel(newdata, ...))
      }
    })
  
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}

# ------------------------------------------------------------------------------
#' add number of cases
#' 
#' Add number of cases in box plot
#'
#' @param y.fun function that calculate position in y-axis
#' @param vjust,position,na.rm refer to \code{\link{stat_summary}}
#'
#' @references \url{https://stackoverflow.com/questions/28846348/add-number-of-observations-per-group-in-ggplot2-boxplot}
#'
#' @examples
#' tibble(x = gl(2, 50), y= rnorm(100)) %>% 
#'   ggplot(aes(x, y)) + 
#'   geom_boxplot() + 
#'   add_n()
#' 
#' @importFrom ggplot2 position_dodge stat_summary
#' @export
# -----------------------------------------------------------------------------#
add_n <- function(
  y.fun= median,
  vjust= -0.3,
  position= position_dodge(width= 0.75),
  na.rm= T
) {
  stat_summary(
    fun.data= function(x) c(y= y.fun(x), label= length(x)),
    geom= 'text',
    vjust= vjust,
    position= position,
    na.rm= na.rm
  )
}
