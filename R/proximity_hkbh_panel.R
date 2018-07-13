#' @name
#' proximity_hkbh_panel
#'
#' @title
#' Proximity a là Hidalgo, Klinger, Barabási, Hausmann for panel
#' data
#'
#' @description
#' The function computes the so called Proximity, as proposed by Hidalgo,
#' Klinger, Barab{\'a}si, and Hausmann in 2007 on \emph{Science}, for a panel
#' data set.
#'
#' @details
#'
#'
#' @references
#' Hidalgo, Klinger, Barab{\'a}si and Hausmann (2007) ``The Product Space
#' Conditions the Development of Nations'', \emph{Science}, 317, 482--487.
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
#' @return A data.frame with the Coherence Index of each geographical area.
#'
#' @examples
#' RCI <- coherence(data = df, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

proximity_hkbh_panel <- function (data, geo_dim, kng_dim, kng_nbr, time_dim,
                                  binary_mode = "RTA") {

    data <- as.data.frame(data)

    time_dim <- deparse(substitute(time_dim))
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    time_span <- unique(data[, time_dim])
    Phi <- lapply(time_span, function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]

        pp <- proximity_hkbh(data_subset, geo_dim, kng_dim, kng_nbr,
                             binary_mode, names_as_strings = TRUE)

        diversification <- attr(pp, 'diversification')
        diversification <- cbind(geo_dim = names(diversification),
                                 time_dim = t,
                                 diversification)
        ubiquity <- attr(pp, 'ubiquity')
        ubiquity <- cbind(kng_dim = names(ubiquity),
                          time_dim = t,
                          ubiquity)

        colnames(pp) <- c(paste0(kng_dim, "_1"),
                          paste0(kng_dim, "_2"),
                          as.character(t))

        return(list(pp,
                    diversification, ubiquity))
    })


    diversification <- sapply(Phi, "[", 2)
    diversification <- do.call("rbind.data.frame", diversification)

    ubiquity <- sapply(Phi, "[", 3)
    ubiquity <- do.call("rbind.data.frame", ubiquity)

    Phi <- sapply(Phi, "[", 1)
    Phi <- plyr::join_all(Phi)
    Phi <- reshape2::melt(Phi)
    colnames(Phi) <- c(paste0(kng_dim, "_1"),
                       paste0(kng_dim, "_2"),
                       time_dim,
                       "Phi")

    class(Phi) <- c('data.frame', 'reks_proximity_hkbh')
    attr(Phi, 'diversification') <- diversification
    attr(Phi, 'ubiquity') <- ubiquity
    if (binary_mode == 'RTA') {
        attr(Phi, "binary_mode") <- 'RTA'
    }
    if (binary_mode == 'RCA') {
        attr(Phi, "binary_mode") <- 'RCA'
    }
    if (binary_mode == 'simple') {
        attr(Phi, "binary_mode") <- 'simple'
    }
    if (binary_mode == "higher_quartiles") {
        attr(Phi, "binary_mode") <- 'higher_quartiles'
    }
    if (binary_mode == "higher_quartiles_kng") {
        attr(Phi, "binary_mode") <- 'higher_quartiles_kng'
    }
    attr(Phi, "geo_dim") <- geo_dim
    attr(Phi, "kng_dim_1") <- paste0(kng_dim, "_1")
    attr(Phi, "kng_dim_2") <- paste0(kng_dim, "_2")
    attr(Phi, "kng_nbr") <- kng_nbr
    attr(Phi, "time_dim") <- time_dim
    attr(Phi, "prx") <- "Phi"
    attr(Phi, "directed") <- FALSE

    return(Phi)
}
