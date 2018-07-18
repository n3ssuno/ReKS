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
#' @param data It is expected to be a dataframe in "long" format.
#' @param geo_dim It is the name of the column of the data.frame that
#' represents its geographical dimension (e.g., the different regions of
#' analysis).
#' @param kng_dim It is the name of the column of the data.frame that
#' represents its knowledge dimension (e.g., the different patent classes of
#' analysis).
#' @param kng_nbr It is the name of the column of the data.frame that
#' represents the numerosity of each knowledge class (e.g., the number of
#' patents a region has in a given year in a given patent class).
#' @param null_model It is useful to reuse this code in the coherence_panel
#' function, leave if NULL unless you have an object of class reks_null_model
#' that you whant to use as null model.
#' @param names_as_strings If TRUE the function assumes you have written the
#' name of the columns above said as strings, otherwise they are assumed to
#' be symbols. Default is FALSE.
#' @return A data.frame with the Coherence Index of each geographical area.
#'
#' @examples
#' RCI <- coherence(data = df, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

coherence <- function(data, geo_dim, kng_dim, kng_nbr,
                      null_model = NULL, names_as_strings = FALSE) {

    # Preliminary operations, checks and transformations ---------------

    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop(paste0("Package \"reshape2\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
        stop(paste0("Package \"plyr\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    data <- as.data.frame(data)

    if (names_as_strings == FALSE) {
        geo_dim <- deparse(substitute(geo_dim))
        kng_dim <- deparse(substitute(kng_dim))
        kng_nbr <- deparse(substitute(kng_nbr))
    }

    data <- data[complete.cases(data[, c(geo_dim, kng_dim, kng_nbr)]), ]

    # Null model (Hypergeometric distribution) ---------------

    if (is.null(null_model)) {
        C <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr, 'simple')

        J <- t(C) %*% C

        n <- colSums(C)

        mu <- n %*% t(n)
        mu <- mu/nrow(C)
        rownames(mu) <- colnames(mu)

        s2 <- (1 - n/nrow(C)) %*% t((nrow(C) - n)/(nrow(C) - 1))
        s2 <- mu * s2
        rownames(s2) <- colnames(s2)

        t <- (J - mu)/sqrt(s2)

        diag(t) <- 0

        null_model <- list()
        null_model$t <- t
        class(null_model) <- "reks_null_model"
    }

    # Coherence index ---------------

    if (class(null_model) != "reks_null_model") {
        stop('null_model must be of "reks_null_model" class')
    }

    # ee <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr, 'simple')
    ee <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr, 'none')

    en <- colnames(ee)
    nc <- setdiff(colnames(null_model$t), colnames(ee))
    en <- c(en, nc)
    for (i in nc) {
        ee <- cbind(ee, 0)
    }
    colnames(ee) <- en
    col.order <- colnames(null_model$t)
    ee <- ee[, col.order]

    ones <- matrix(1, nrow = nrow(null_model$t), ncol = nrow(null_model$t))
    diag(ones) <- 0

    WAR <- ee %*% null_model$t / ee %*% ones
    R <- rowSums(WAR * ee / rowSums(ee))
    R[is.nan(R)] <- 0
    R <- as.data.frame(R)
    R <- cbind(rownames(R), R)
    measure <- "Coherence"
    colnames(R) <- c(geo_dim, measure)

    class(R) <- c("reks_coherence", "data.frame")
    attr(R, 'geo_dim') <- geo_dim
    attr(R, 'kng_dim') <- kng_dim
    attr(RCI, 'measure') <- measure
    attr(R, 'null_model') <- null_model

    return(R)
}
