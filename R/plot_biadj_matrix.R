#' @name
#' plot_biadj_matrix
#'
#' @title
#' Plot a biadjacency matrix
#'
#' @description
#' The function plots a biadjacency matrix.
#'
#' @details
#' The function transforms the ``long'' table provided in a biadjacency matrix
#' and than plot it.
#'
#' @encoding UTF-8
#'
#' @param data It must be a ``long'' dataframe.
#' @param geo_dim It must be a symbolic variable that denotes the name of the
#' column of the dataframe that represents the geographical units of analysis.
#' @param kng_dim It must be a symbolic variable that denotes the name of the
#' column of the dataframe that represents the technological knowledge units of
#' analysis.
#' @param kng_nbr It must be a symbolic variable that denotes the name of the
#' column of the dataframe that represents the number of technological units
#' of a given type that a given geographical units owns.
#' @param binary_mode It must be a string that represents the way in which the
#' data are transformed from count values to binary values. It can be 'none',
#' 'simple', 'RTA', 'RCA', 'higher_quartiles', 'higher_quartiles_kng',
#' 'higher_deciles_kng'.
#' @param order It must be a string that represents the way in which the matrix
#' rows and columns are re-arranged. By now the only possible option is "DU" that arranges the rows
#' and columns in order inverse to their degree (i.e., to the
#' \emph{diversification} of the geographical area and the \emph{ubiquity} of
#' the knowledge dimension considered.)
#' @param ... Other parameters taken by plot.
#' @return Nothing

plot_biadj_matrix <- function(data, geo_dim, kng_dim, kng_nbr,
                              binary_mode, order = "DU", ...) {
    # Order can be
    # - "DU" = diversification - ubiquity
    # - "FC" = fitness - complexity
    # - "HH" = Hidalgo-Hausmann complexity
    # TODO
    # The last two haven't been implemented yet
    data <- as.data.frame(data)
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))
    BM <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr, binary_mode)
    plot(BM, order, ...)
}


# https://www.r-bloggers.com/creating-an-image-of-a-matrix-in-r-using-image/

#' @name
#' plot.rks_biadj_matrix
#'
#' @title
#' Plot a biadjacency matrix
#'
#' @description
#' The function plots a biadjacency matrix.
#'
#' @details
#' Use \code{\link{plot_biadj_matrix}} instead of this function to have a better
#' and easier user experience.
#'
#' @encoding UTF-8
#'
#' @param x It is expected to be an object of class \code{rks_biadj_matrix}.
#' @param order It is the way in which rows and columns of the matrix are
#' re-arranged. By now the only possible option is "DU" that arranges the rows
#' and columns in order inverse to their degree (i.e., to the
#' \emph{diversification} of the geographical area and the \emph{ubiquity} of
#' the knowledge dimension considered.)
#' @param ... Other parameters taken by plot.
#' @return Nothing
#'
#' @export

plot.rks_biadj_matrix <- function(x, order = "DU", ...) {
    if (order == "DU") {
        du <- .get_du(x)
        row_order <- order(du$diversification)
        col_order <- order(du$ubiquity, decreasing = T)
    }
    x <- x[row_order, col_order]
    rotate <- function(m) t(apply(m, 2, rev))
    image(rotate(x),
          col = c('white', 'black'),
          xlab = attr(x, "kng_dim"), ylab = attr(x, "geo_dim"),
          axes = FALSE, ...)
    box()
    invisible()
}
