#' @name
#' rta
#'
#' @title
#' Relative Technological Advantages
#'
#' @description
#' Returns a matrix with the Relative Technological Advantages values for each
#' geo. area for each tech. class (or a binary matrix, if specified).
#'
#' @details
#' Returns a matrix with the Relative Technological Advantages values for each
#' geo. area for each tech. class. If specified, the matrix is made binary
#' putting a 1 (TRUE) each time the RTAs are above 1, and 0 otherwise.
#'
#' @encoding UTF-8
#'
#' @param data It is expected to be a matrix (use occurrence_matrix or
#' cooccurrence_matrix to get it from a "long" data.frame)
#' @param binary
#' @return A sparse matrix

rta <- function(data, binary = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop(paste0("Package \"Matrix\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    RA <- Matrix::t(
        Matrix::t(
            data / Matrix::rowSums(data)) /
            (Matrix::colSums(data) / sum(data)))

    if (isTRUE(binary)) {
        RA <- as(RA, "ngCMatrix")
    }

    return(RA)
}
