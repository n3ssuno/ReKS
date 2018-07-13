#' @name
#' proximity_zctp_panel
#'
#' @title
#' Proximity a là Zaccaria, Cristelli, Tacchella, Pietronero, for
#' panel data
#'
#' @description
#' The function computes the so called Proximity, as proposed by Zaccaria,
#' Cristelli, Tacchella, and Pietronero in 2014 on \emph{PLoS ONE}, for a panel
#' data set.
#'
#' @details
#'
#'
#' @references
#' Zaccaria, Cristelli, Tacchella and Pietronero (2014) ``How the Taxonomy of
#' Products Drives the Economic Development of Countries'', \emph{PLoS ONE},
#' 0113770, 1--17.
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

proximity_zctp_panel <- function (data, geo_dim, kng_dim, kng_nbr, time_dim,
                                  binary_mode = "RTA") {

    data <- as.data.frame(data)

    time_dim <- deparse(substitute(time_dim))
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    time_span <- unique(data[, time_dim])
    B <- lapply(time_span, function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]

        pp <- proximity_zctp(data_subset, geo_dim, kng_dim, kng_nbr,
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

    diversification <- sapply(B, "[", 2)
    diversification <- do.call("rbind.data.frame", diversification)

    ubiquity <- sapply(B, "[", 3)
    ubiquity <- do.call("rbind.data.frame", ubiquity)

    B <- sapply(B, "[", 1)
    B <- plyr::join_all(B)
    B <- reshape2::melt(B)
    colnames(B) <- c(paste0(kng_dim, "_1"),
                       paste0(kng_dim, "_2"),
                       time_dim,
                       "B")

    class(B) <- c('data.frame', 'reks_proximity_zctp')
    attr(B, 'diversification') <- diversification
    attr(B, 'ubiquity') <- ubiquity
    if (binary_mode == 'RTA') {
        attr(B, "binary_mode") <- 'RTA'
    }
    if (binary_mode == 'RCA') {
        attr(B, "binary_mode") <- 'RCA'
    }
    if (binary_mode == 'simple') {
        attr(B, "binary_mode") <- 'simple'
    }
    if (binary_mode == "higher_quartiles") {
        attr(B, "binary_mode") <- 'higher_quartiles'
    }
    if (binary_mode == "higher_quartiles_kng") {
        attr(B, "binary_mode") <- 'higher_quartiles_kng'
    }
    attr(B, "geo_dim") <- geo_dim
    attr(B, "kng_dim_1") <- paste0(kng_dim, "_1")
    attr(B, "kng_dim_2") <- paste0(kng_dim, "_2")
    attr(B, "kng_nbr") <- kng_nbr
    attr(B, "time_dim") <- time_dim
    attr(B, "prx") <- "B"
    attr(B, "directed") <- TRUE

    return(B)
}
