#' @name
#' coherence
#'
#' @title
#' Regional Coherence Index
#'
#' @description
#' The function computes the so called Coherence index.
#'
#' @details
#' The function computes the so called Coherence index.
#' It assumes that the "universe" from which you derive the distribution of
#' reference is composed by all the information provided in the database, and
#' only this.
#' You can use it both on a panel data set (if you identify also a column with
#' a temporal indication of the observations) or on cross-section data (by
#' leaveing the parameter just said as NULL).
#'
#' @references
#' Teece, Rumelt, Dosi and Winter (1994) ``Understanding Corporate
#' Coherence: Theory and Evidence'', \emph{Journal of Economic Behavior &
#' Organization}, 23, 1-30;
#'
#' Nesta and Saviotti (2005) ``Coherence of the Knowledge Base and the
#' Firm's Innovative Performance: Evidence from the U.S. Pharmaceutical
#' Industry'', \emph{Journal of Industrial Economics}, 53, 123-142;
#'
#' Nesta and Saviotti (2006) ``Firm Knowledge and Market Value in
#' Biotechnology'', \emph{Industrial and Corporate Change}, 15, 625-652;
#'
#' Bottazzi and Pirino (2010) ``Measuring Industry Relatedness and
#' Corporate Coherence'', \emph{SSRN Electronic Journal}, 11, 1--24;
#'
#' Quatraro (2010) ``Knowledge Coherence, Variety and Economic Growth:
#' Manufacturing Evidence from Italian Regions'', \emph{Research Policy}, 39,
#' 1289-1302;
#'
#' Rocchetta and Mina (2017), ``Technological Coherence and the Adaptive
#' Resilience of Regional Economies'', \emph{Papers in Evolutionary Economic
#' Geography}, Utrecht University.
#'
#' @encoding UTF-8
#'
#' @param data It is expected to be a matrix (use occurrence_matrix to get it
#' from a "long" data.frame)
#' @param relatedness
#' @return A data.frame with the Coherence Index of each geographical area.
#'
#' @examples
#' RCI <- coherence(data = df, relatedness = relMtx)

coherence <- function(occurrence_mtx, relatedness_mtx) {
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop(paste0('Package \"Matrix\" needed for this function to work. ',
                    'Please install it.'), call. = FALSE)

    geo_dim <- attr(occurrence_mtx, 'geo_dim')
    kng_dim <- attr(occurrence_mtx, 'kng_dim')
    if (is.list(occurrence_mtx))
        time_dim <- attr(occurrence_mtx, 'time_dim')
    measure <- "Coherence"

    coherence_crossSection <- function(occurrence_mtx, relatedness_mtx) {

        # Preliminary operations, checks and transformations ---------------

        oc_mtx_names <- colnames(occurrence_mtx)
        rl_mtx_names <- colnames(relatedness_mtx)
        if (dim(occurrence_mtx)[[2]] != sum(dim(relatedness_mtx)) / 2) {
            names_tbr <- setdiff(rl_mtx_names, oc_mtx_names)
            if (length(names_tbr) != 0)
                relatedness_mtx <- relatedness_mtx[
                    -which(rownames(relatedness_mtx) %in% names_tbr),
                    -which(colnames(relatedness_mtx) %in% names_tbr)]
        }
        if (dim(occurrence_mtx)[[2]] != sum(dim(relatedness_mtx)) / 2)
            stop(paste('There is some problem, because the two matrices',
                       'considered have a different number of',
                       'columns.'), call. = FALSE)

        rl_mtx_names <- colnames(relatedness_mtx)
        if (any(oc_mtx_names != rl_mtx_names)) {
            oc_mtx_names <- oc_mtx_names[, order(colnames(oc_mtx_names))]
            relatedness_mtx <- relatedness_mtx[order(rownames(relatedness_mtx)),
                                               order(colnames(relatedness_mtx))]
        }
        if (any(oc_mtx_names != rl_mtx_names))
            stop(paste('There is some problem, because the there is no perfect',
                       'correspondence between the column names of the two',
                       'matrices considered.'), call. = FALSE)

        ones <- !Matrix::diag(TRUE,
                              nrow = nrow(relatedness_mtx),
                              ncol = nrow(relatedness_mtx))

        # Waighted Average Relatedness

        WAR_num <- Matrix::tcrossprod(occurrence_mtx, relatedness_mtx)
        WAR_den <- Matrix::tcrossprod(occurrence_mtx, ones)
        WAR <- WAR_num / WAR_den

        # Coherence

        C_num <- WAR * occurrence_mtx
        C_den <- Matrix::rowSums(occurrence_mtx)
        C <- Matrix::rowSums(C_num / C_den)
        C[which(is.nan(C))] <- 0

        C <- cbind.data.frame(names(C), unlist(C))
        colnames(C) <- c(geo_dim, measure)

        return(C)
    }

    coherence_panel <- function(occurrence_mtx, relatedness_mtx) {

        time_span <- names(occurrence_mtx)
        C <- lapply(as.character(time_span),
                    function(y) coherence_crossSection(occurrence_mtx[[y]],
                                                       relatedness_mtx))
        yrs <- unlist(mapply(rep, time_span, lapply(C, nrow)))
        C <- do.call("rbind", C)
        C <- cbind.data.frame(C, yrs)
        C <- C[, c(1, 3, 2)]
        colnames(C) <- c(geo_dim, time_dim, measure)

        return(C)
    }

    fntn <- ifelse(is.list(occurrence_mtx),
                   "coherence_panel",
                   "coherence_crossSection")
    C <- do.call(fntn, list(occurrence_mtx, relatedness_mtx))

    # final steps

    # class(R) <- c("reks_coherence", "data.frame")
    attr(C, 'geo_dim') <- geo_dim
    attr(C, 'kng_dim') <- kng_dim
    if (is.list(occurrence_mtx))
        attr(C, 'time_dim') <- time_dim
    attr(C, 'measure') <- measure

    return(C)
}

