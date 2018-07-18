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
#' @param binary_mode It is "RTA" by default. The possible values are: "simple"
#' (if there is at least a patent in the geo. area in the particular tech. class
#' it is 1, and 0 otherwise); "RTA" or "RCA" (Balassa method);
#' "higher_quartiles" (similar to "simple", but it exludes the first quartile
#' (it seems the best choice, see below); "higher_quartiles_kng" (similar to
#' 'higher_quartiles', but the quartiles are computed for each knowledge class,
#' assuming that some are more "ubiquitous" than others); "higher_deciles_kng"
#' (similar to 'higher_quartiles_kng', but it excludes the lower 10% rather than
#' the lower 25%). Warning: look carefully at the plot_biadj_matrix() output
#' because I personally observed that with RTA applied to patent data you do not
#' have the triangular structure observed in the trade data.
#' @param scale It is TRUE by default. Otherwise, the Complexity Index is not
#' standardised (CI - mean[CI] / sd[CI]).
#' @param names_as_strings If TRUE the function assumes you have written the
#' name of the columns above said as strings, otherwise they are assumed to
#' be symbols. Default is FALSE.
#' @return A data.frame with the Complexity Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RKCI <- complexity(data = df, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

complexity_hh <- function(data, geo_dim, kng_dim, kng_nbr,
                          binary_mode = "RTA", scale = TRUE,
                          names_as_strings = FALSE) {
    # if (!requireNamespace("reshape2", quietly = TRUE)) {
    #     stop(paste0("Package \"reshape2\" needed for this function to work. ",
    #                 "Please install it."), call. = FALSE)
    # }

    data <- as.data.frame(data)

    if (names_as_strings == FALSE) {
        geo_dim <- deparse(substitute(geo_dim))
        kng_dim <- deparse(substitute(kng_dim))
        kng_nbr <- deparse(substitute(kng_nbr))
    }

    data <- data[complete.cases(data[, c(geo_dim, kng_dim, kng_nbr)]), ]

    mm <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr,
                                binary_mode)

    if (any(rowSums(mm) == 0)) {
        mm <- mm[-which(rowSums(mm) == 0), ]
    }
    if (any(colSums(mm) == 0)) {
        mm <- mm[, -which(colSums(mm) == 0)]
    }

    du <- .get_du(mm)

    mm_tilde <- t(t(mm) / du$ubiquity) %*% t(mm)
    mm_tilde <- mm_tilde / du$diversification

    if (!all(round(rowSums(mm_tilde)) == 1)) {
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

    RKCI <- cbind.data.frame(rownames(mm_tilde), RKCI)
    measure <- "Complexity"
    colnames(RKCI) <- c(geo_dim, measure)

    if (scale == TRUE) {
        RKCI[, measure] <- scale(RKCI[, measure])
        warning('The values of the index have been standardised.')
    }


    class(RKCI) <- c('data.frame', 'reks_hh_complexity')
    attr(RKCI, 'diversification') <- du$diversification
    attr(RKCI, 'ubiquity') <- du$ubiquity
    attr(RKCI, 'standardised') <- ifelse(scale == TRUE, TRUE, FALSE)
    if (binary_mode == 'RTA') {
        attr(RKCI, "binary_mode") <- 'RTA'
    }
    if (binary_mode == 'RCA') {
        attr(RKCI, "binary_mode") <- 'RCA'
    }
    if (binary_mode == 'simple') {
        attr(RKCI, "binary_mode") <- 'simple'
    }
    if (binary_mode == "higher_quartiles") {
        attr(RKCI, "binary_mode") <- 'higher_quartiles'
    }
    if (binary_mode == "higher_quartiles_kng") {
        attr(RKCI, "binary_mode") <- 'higher_quartiles_kng'
    }
    attr(RKCI, 'geo_dim') <- geo_dim
    attr(RKCI, 'kng_dim') <- kng_dim
    attr(RKCI, 'measure') <- measure

    return(RKCI)
}
