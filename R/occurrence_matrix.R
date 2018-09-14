#' @name
#' occurrence_matrix
#'
#' @title
#' Occurrences matrix
#'
#' @description
#' Returns an occurrences matrix.
#'
#' @details
#' The function returns the occurrences matrix (or a list of matrices,
#' in case of panel data), based on a data.frame in "long" format.
#'
#' Possible binary_mode options:
#' \itemize{
#'  \item{"none"}{weighted matrix (not binary);}
#'  \item{"simple"}{if there is at least a patent in the geo. area
#'    in the particular tech. class it is 1, and 0 otherwise;}
#'  \item{"RTA"}{Balassa method;}
#'  \item{"RCA"}{Balassa method;}
#' }
#'
#' @encoding UTF-8
#'
#' @param data It is expected to be a data.frame in "long" format.
#' @param geo_dim It is the name of the column of the data.frame that
#' represents its geographical dimension (e.g., the different regions of
#' analysis).
#' @param kng_dim It is the name of the column of the data.frame that
#' represents its knowledge dimension (e.g., the different patent classes of
#' analysis).
#' @param kng_nbr It is the name of the column of the data.frame that
#' represents the numerosity of each knowledge class (e.g., the number of
#' patents a region has in a given year in a given patent class).
#' @param time_dim
#' @param binary_mode
#' @return A sparse matrix (or a list of matrices).
#'
#' @examples
#' ocMtx <- occurrence_matrix(data = df,
#'            geo_dim = NUTS2, kng_dim = IPC3,
#'            kng_nbr = Npatents)
#'
#' ocMtx <- occurrence_matrix(data = df,
#'            geo_dim = NUTS2, kng_dim = IPC3,
#'            kng_nbr = Npatents,
#'            time_dim = "Year")
#'
#' ocMtx <- occurrence_matrix(data = df,
#'            geo_dim = NUTS2, kng_dim = IPC3,
#'            binary_mode = "simple")

occurrence_matrix <- function(data, geo_dim, kng_dim,
                              kng_nbr = NULL,
                              time_dim = NULL,
                              binary_mode = "none") {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop(paste0("Package \"Matrix\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    # Preliminary controls ---------------

    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))
    time_dim <- deparse(substitute(time_dim))
    if (kng_nbr == "NULL")
        kng_nbr <- NULL
    if (time_dim == "NULL")
        time_dim <- NULL

    if (binary_mode != 'simple' & is.null(kng_nbr)) {
        stop('Either you specify that you want a \"simple\" matrix, or you have
             to specify a column for the number of pieces of knowledge')
    }

    data <- as.data.frame(data)

    if (is.null(kng_nbr)) {
        data <- unique(data[, c(geo_dim, kng_dim, time_dim)])
    } else {
        if (is.null(time_dim)) {
            if (anyDuplicated(data[, c(geo_dim, kng_dim)])) {
                frml <- formula(paste(kng_nbr, "~",
                                      geo_dim, "+", kng_dim))
                data <- aggregate(formula = frml, data = data, FUN = sum)
                warning(paste('Since there are duplicated cases,',
                              'the function has collapsed them,',
                              'by summing the number of pieces of knowledge'),
                        call. = F)
            }
        } else {
            if (anyDuplicated(data[, c(geo_dim, kng_dim, time_dim)])) {
                frml <- formula(paste(kng_nbr, "~",
                                      geo_dim, "+", kng_dim, "+", time_dim))
                data <- aggregate(formula = frml, data = data, FUN = sum)
                warning(paste('Since there are duplicated cases,',
                              'the function has collapsed them,',
                              'by summing the number of pieces of knowledge'),
                        call. = F)
            }
        }
    }

    get_mtx <- function(data, frml, binary_mode) {
        bm <- xtabs(formula = formula(frml),
                    data = data,
                    sparse = TRUE)

        bm <- switch(binary_mode,
                     RTA = rta(bm, binary = TRUE),
                     RCA = rta(bm, binary = TRUE),
                     simple = as(bm, "ngCMatrix"),
                     none = bm)

        return(bm)
    }

    # Main function ---------------

    frml <- ifelse(is.null(kng_nbr),
                   paste("~", geo_dim, "+", kng_dim),
                   paste(kng_nbr, "~", geo_dim, "+", kng_dim))
    if (is.null(time_dim))
        BM <- get_mtx(data, frml, binary_mode)
    else {
        time_span <- unique(data[, time_dim])
        BM <- lapply(time_span, function(y)
            get_mtx(data[which(data[, time_dim] == y), ], frml, binary_mode))
        names(BM) <- time_span
    }

    # Closing operations

    # class(BM) <- c('reks_biadj_matrix', 'matrix')
    attr(BM, "geo_dim") <- geo_dim
    attr(BM, "kng_dim") <- kng_dim
    if (!is.null(time_dim))
        attr(BM, "time_dim") <- time_dim
    attr(BM, "binary")  <- ifelse(binary_mode == "none", FALSE, TRUE)
    if (binary_mode == 'RTA')
        attr(BM, "binary_mode") <- 'RTA'
    if (binary_mode == 'RCA')
        attr(BM, "binary_mode") <- 'RCA'
    if (binary_mode == 'simple')
        attr(BM, "binary_mode") <- 'simple'

    return(BM)
}
