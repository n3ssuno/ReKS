#' @title equal_dim
#' @param mtx_smaller Matrix with less rows or columns
#' @param mtx_bigger Matrix with more rows or columns
#' @param dim The dimension on which must be one of 'row' (i.e., 1), 'col'
#' (i.e., 2), or 'both' (i.e., \code{c(1, 2)}).
#' @return The smaller matrix with the same dimension of the bigger one. The
#' additional rows or columns are zeros.
#' @keywords internal

.equal_dim <- function(mtx_smaller, mtx_bigger, MARGIN = NULL) {
    # Add rows or columns to mtx_smaller (zeros) so that
    #  it has the same rows and columns of mtx_bigger
    if (is.null(MARGIN)) {
        stop(paste0("MARGIN must be one of 'row', 'col', or 'both'.
                    And you can also use the 'apply' notation."))
    }
    if (any(MARGIN %in% c('row', 1, 'both'))) {
        rns2 <- setdiff(rownames(mtx_bigger), rownames(mtx_smaller))
        nms <- rownames(mtx_smaller)
        mtx_smaller <- rbind(mtx_smaller, matrix(0,
                                                 nrow = length(rns2),
                                                 ncol = ncol(mtx_smaller)))
        rownames(mtx_smaller) <- c(nms, rns2)
        mtx_smaller <- mtx_smaller[match(rownames(mtx_bigger),
                                         rownames(mtx_smaller)), ]
    }
    if (any(MARGIN %in% c('col', 2, 'both'))) {
        cns2 <- setdiff(colnames(mtx_bigger), colnames(mtx_smaller))
        nms <- colnames(mtx_smaller)
        mtx_smaller <- cbind(mtx_smaller, matrix(0,
                                                 nrow = nrow(mtx_smaller),
                                                 ncol = length(cns2)))
        colnames(mtx_smaller) <- c(nms, cns2)
        mtx_smaller <- mtx_smaller[, match(colnames(mtx_bigger),
                                           colnames(mtx_smaller))]
    }

    return(mtx_smaller)
}

#' @name .get_ji
#' @param mtx1 First matrix
#' @param mtx2 Second matrix
#' @return Jaccard similarity coefficient
#' @keywords internal

.get_ji <- function(mtx1, mtx2) {
    if (nrow(mtx1) < nrow(mtx2)) {
        mtx1 <- .equal_dim(mtx1, mtx2, 'row')
    }
    if (nrow(mtx1) > nrow(mtx2)) {
        mtx2 <- .equal_dim(mtx2, mtx1, 'row')
    }
    if (ncol(mtx1) < ncol(mtx2)) {
        mtx1 <- .equal_dim(mtx1, mtx2, 'col')
    }
    if (ncol(mtx1) > ncol(mtx2)) {
        mtx2 <- .equal_dim(mtx2, mtx1, 'col')
    }
    tb <- sum(mtx1 == 1 & mtx2 == 1)
    t1 <- sum(mtx1 == 1)
    t2 <- sum(mtx2 == 1)
    JI <- tb / (t1 + t2 - tb)

    return(JI)
}

#' @name jaccard_index
#'
#' @title
#' Jaccard similarity coefficient
#'
#' @description
#' The function computes the Jaccard similarity coefficient.
#'
#' @details
#' The function computes the Jaccard similarity coefficient.
#'
#' @references
#'
#' @encoding UTF-8
#'

jaccard_index <- function(data, geo_dim, kng_dim, kng_nbr, time_dim,
                          binary_mode = "RTA") {

    data <- as.data.frame(data)

    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))
    time_dim <- deparse(substitute(time_dim))

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    time_span <- unique(data[, time_dim])
    mtx <- lapply(time_span, function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]
        m <- .get_biadj_matrix(data_subset,
                               geo_dim, kng_dim, kng_nbr,
                               binary_mode)
        return(m)
    })

    JI <- lapply(1:(length(time_span) - 1), function(t) {
        c(firstYear = time_span[[t]],
          sectonYear = time_span[[t + 1]],
          JI = .get_ji(mtx[[t]], mtx[[t + 1]]))
    })
    JI <- do.call('rbind', JI)
    JI <- as.data.frame(JI)

    class(JI)   <- c('reks_jaccard_index', 'data.frame')

    return(JI)
}

#' @name plot.reks_jaccard_index
#' @param JI Object of class \bold{plot.reks_jaccard_index}
#' @export

plot.reks_jaccard_index <- function(JI) {
    plot(JI[, "firstYear"], JI[, "JI"], type = "b",
         xlab = 'Year', ylab = 'Jaccard index')
}
