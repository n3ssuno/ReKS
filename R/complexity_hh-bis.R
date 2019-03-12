#' @name
#' complexity_hh
#'
#' @title
#' Regional Knowledge Complexity Index a là Hidalgo-Hausmann
#'
#' @description
#' The function computes the ``Technological Complexity'' of each geographical
#' area considered, in each given year, in line with the so called \emph{Method
#' of Reflections} (see the references below).
#'
#' @references
#' Hidalgo and Hausmann (2009) ``The Building Blocks of Economic
#' Complexity'', \emph{PNAS}, 106, 10570--10575;
#'
#' Antonelli, Crespi, Mongeau Ospina and Scellato (2017) ``Knowledge
#' Composition, Jacobs Externalities and Innovation Performance in European
#' Regions'', \emph{Regional Studies}, 51, 1708--1720;
#'
#' Balland and Rigby (2017) ``The Geography of Complex Knowledge'',
#' \emph{Economic Geography}, 93, 1--23.
#'
#' @encoding UTF-8
#'
#' @param occurrence_mtx It is expected to be an object created by the
#' ReKS::occurrence_matrix function.
#' @param rta If TRUE (default) it uses the RTA/RCA of the original matrix.
#' @param binary If TRUE (default) it binarises the RTA/RCA matrix
#' (can be used only together with rta=TRUE).
#' @param scale If TRUE (default) the output is standardised
#' (z-scores = CI - mean[CI] / sd[CI]).
#' @return A data.frame with the Complexity Index of each geographical area
#' in each time step of analysis.
#' @examples
#' OccM <- occurrence_matrix(data = DATA,
#'            geo_dim = NUTS2, kng_dim = IPC4, kng_nbr = N_Patents)
#' RKCI <- complexity(OccM)

complexity_hh <- function(occurrence_mtx,
                          rta = TRUE, binary = TRUE, scale = TRUE) {
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop(paste0('Package \"Matrix\" needed for this function to work. ',
                    'Please install it.'), call. = FALSE)

    geo_dim <- attr(occurrence_mtx, 'geo_dim')
    kng_dim <- attr(occurrence_mtx, 'kng_dim')
    if (is.list(occurrence_mtx))
        time_dim <- attr(occurrence_mtx, 'time_dim')
    measure <- "Complexity"

    complexity_hh_crossSection <- function(occurrence_mtx) {
        if (any(Matrix::rowSums(occurrence_mtx) == 0)) {
            occurrence_mtx <- occurrence_mtx[-Matrix::which(
                Matrix::rowSums(occurrence_mtx) == 0), ]
        }
        if (any(Matrix::colSums(occurrence_mtx) == 0)) {
            occurrence_mtx <- occurrence_mtx[, -Matrix::which(
                Matrix::colSums(occurrence_mtx) == 0)]
        }

        rnms <- rownames(occurrence_mtx)
        # cnms <- colnames(occurrence_mtx)

        if (isTRUE(rta)) {
            occurrence_mtx <- rta(occurrence_mtx)
            if (isTRUE(binary))
                occurrence_mtx <- Matrix::Matrix(
                    ifelse(occurrence_mtx >= 1, 1, 0),
                    nrow = nrow(occurrence_mtx))
        }

        du <- ReKS:::.get_du(Matrix::as.matrix(occurrence_mtx))
        #TODO

        mm_tilde <- Matrix::t(Matrix::t(occurrence_mtx) / du$ubiquity)
        mm_tilde <- Matrix::tcrossprod(mm_tilde, occurrence_mtx)
        mm_tilde <- mm_tilde / du$diversification

        if (!all(round(Matrix::rowSums(mm_tilde)) == 1)) {
            stop(paste("The matrix is not row-stochastic.\n",
                       "It is not possible to compute the measure."))
        }
        if (round(Re(as.complex(eigen(mm_tilde)$value[1]))) != 1) {
            stop(paste("The first eigen-value is different from 1.",
                       "It is not possible to compute the measure."))
        }

        RKCI <- eigen(mm_tilde)$vectors
        if (dim(RKCI)[2] >= 2) {
            RKCI <- RKCI[, 2]
            RKCI <- Re(as.complex(RKCI))

            if (cor(RKCI, du$diversification,
                    use = "pairwise.complete.obs", method = "spearman") < 0) {
                RKCI <- -RKCI
            }
        } else {
            RKCI <- NA
        }

        RKCI <- cbind.data.frame(rnms,
                                 RKCI)
        colnames(RKCI) <- c(geo_dim, measure)

        if (scale == TRUE) {
            RKCI[, measure] <- scale(as.numeric(RKCI[, measure]))
            # warning('The values of the index have been standardised.')
        }

        return(RKCI)
    }

    complexity_hh_panel <- function(occurrence_mtx) {
        time_span <- names(occurrence_mtx)
        RKCI <- lapply(as.character(time_span),
                    function(y) complexity_hh_crossSection(occurrence_mtx[[y]]))
        yrs <- unlist(mapply(rep, time_span, lapply(RKCI, nrow)))
        RKCI <- do.call("rbind", RKCI)
        RKCI <- cbind.data.frame(RKCI, yrs)
        RKCI <- RKCI[, c(1, 3, 2)]
        colnames(RKCI) <- c(geo_dim, time_dim, measure)

        return(RKCI)
    }

    fntn <- ifelse(is.list(occurrence_mtx),
                   "complexity_hh_panel",
                   "complexity_hh_crossSection")
    Cx <- do.call(fntn, list(occurrence_mtx))

    # final steps

    # class(R) <- c("reks_complexity_hh", "data.frame")
    attr(Cx, 'geo_dim') <- geo_dim
    attr(Cx, 'kng_dim') <- kng_dim
    if (is.list(occurrence_mtx))
        attr(Cx, 'time_dim') <- time_dim
    attr(Cx, 'measure') <- measure
    # attr(Cx, 'diversity') <- du$diversification
    # attr(Cx, 'ubiquity') <- du$ubiquity
    attr(Cx, 'standardised') <- scale
    attr(Cx, "RTA") <- rta
    attr(Cx, "binary") <- binary

    return(Cx)
}
