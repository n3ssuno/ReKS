#' @name
#' complexity
#'
#' @title
#' Regional Knowledge Complexity Index
#'
#' @description
#'
#' @details
#'
#' See:
#' \itemize{
#' \item{Hidalgo, Klinger, Barab{\'a}si and Hausmann (2007) "The Product Space
#' Conditions the Development of Nations", \emph{Science}, 317, 482--487;}
#' \item{Hidalgo and Hausmann (2009) "The Building Blocks of Economic
#' Complexity", \emph{PNAS}, 106, 10570--10575;}
#' \item{Antonelli, Crespi, Mongeau Ospina and Scellato (2017) "Knowledge
#' Composition, Jacobs Externalities and Innovation Performance in European
#' Regions", \emph{Regional Studies}, 51, 1708--1720;}
#' \item{Balland and Rigby (2017) "The Geography of Complex Knowledge",
#' \emph{Economic Geography}, 93, 1--23.}
#' }
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
#' @param time_dim It is the name of the column of the data.frame that
#' represents its temporal dimension (e.g., the different years of analysis).
#' It is an optional parameter.
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
#' @return A data.frame with the Complexity Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RKCI <- complexity(data = df, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

complexity <- function(data, geo_dim, kng_dim, kng_nbr, time_dim = NULL,
                       binary_mode = "RTA", scale = TRUE) {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop(paste0("Package \"reshape2\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    data <- as.data.frame(data)

    time_dim <- deparse(substitute(time_dim))
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))

    time_dim_added <- FALSE
    if (time_dim == "NULL") {
        data$AddedTimeDim <- 1
        time_dim <- "AddedTimeDim"
        time_dim_added <- TRUE
    }

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    RKCI <- lapply(unique(data[, time_dim]), function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]
        mm <- .get_biadj_matrix(data_subset, geo_dim, kng_dim, kng_nbr,
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

        CI <- eigen(mm_tilde)$vectors[, 2]
        CI <- Re(as.complex(CI))

        if (cor(CI, du$diversification,
                use = "pairwise.complete.obs", method = "spearman") < 0) {
            CI <- -CI
        }

        CI <- cbind.data.frame(rownames(mm_tilde), t, CI)
        colnames(CI) <- c(geo_dim, time_dim, "Complexity")

        if (scale == TRUE) {
            CI[, 2] <- scale(CI[, 2])
            warning('The values of the index have been standardised.')
        }

        return(CI)
    })

    #names(RKCI) <- unique(data[, time_dim])
    #RKCI <- plyr::join_all(RKCI, geo_dim)
    #RKCI <- reshape2::melt(RKCI)
    RKFI <- do.call("rbind.data.frame", RKFI)
    #colnames(RKCI) <- c(geo_dim, time_dim, "Complexity")
    if (time_dim_added) {
        RKCI <- RKCI[, c(geo_dim, "Complexity")]
    }

    class(RKCI) <- c('data.frame', 'rks_hh_complexity')
    # TODO
    # In that way it returns only the last year
    # attr(RKCI, 'diversification') <- du$diversification
    # attr(RKCI, 'ubiquity') <- du$ubiquity
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

    return(RKCI)
}
