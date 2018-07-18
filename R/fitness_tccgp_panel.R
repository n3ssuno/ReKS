#' @name
#' fitness_tccgp_panel
#'
#' @title
#' Regional Knowledge Fitness Index a là Tacchella, Cristelli,
#' Caldarelli, Gabrielli and Pietronero, for panel data
#'
#' @description
#' The function computes the Fitness (i.e. \emph{competitiveness}) of each given
#' geographical area considered, in each year provided (see the references
#' below for further details).
#'
#' @references
#' Tacchella, Cristelli, Caldarelli, Gabrielli and Pietronero (2012)
#' ``A New Metrics for Countries' Fitness and Products' Complexity'',
#' \emph{Scientific Reports}, 2, 1--7;
#'
#' Tacchella, Cristelli, Caldarelli, Gabrielli and Pietronero (2013)
#' ``Economic Complexity: Conceptual Grounding of a New Metrics for Global
#' Competitiveness'', \emph{Journal of Economic Dynamics and Control}, 37,
#' 1683--1691;
#'
#' Cristelli, Gabrielli, Tacchella, Caldarelli and Pietronero (2013)
#' ``Measuring the Intangibles: A Metrics for the Economic Complexity of
#' Countries and Products'', \emph{PLoS ONE}, 8, e70726;
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
#' @return A data.frame with the Fitness Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RKFI <- fitness(data = df, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

fitness_tccgp_panel <- function(data, geo_dim, kng_dim, kng_nbr, time_dim,
                    binary_mode = "RTA") {
    # if (!requireNamespace("reshape2", quietly = TRUE)) {
    #     stop(paste0("Package \"reshape2\" needed for this function to work. ",
    #                 "Please install it."), call. = FALSE)
    # }

    data <- as.data.frame(data)

    time_dim <- deparse(substitute(time_dim))
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    measure <- "Fitness"

    time_span <- unique(data[, time_dim])
    RKFI <- lapply(time_span, function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]

        FI <- fitness_tccgp(data_subset, geo_dim, kng_dim, kng_nbr,
                            binary_mode, names_as_strings = TRUE)

        iterations <- attr(FI, "iterations")
        convergence <- attr(FI, "convergence")

        diversification <- attr(FI, 'diversification')
        diversification <- cbind(geo_dim = names(diversification),
                                 time_dim = t,
                                 diversification)
        ubiquity <- attr(FI, 'ubiquity')
        ubiquity <- cbind(kng_dim = names(ubiquity),
                          time_dim = t,
                          ubiquity)

        FI <- cbind.data.frame(FI, t)
        FI <- FI[, c(1, 3, 2)]
        colnames(FI) <- c(geo_dim, time_dim, measure)

        return(list(FI,
                    iterations, convergence,
                    diversification, ubiquity))
    })

    iterations <- sapply(RKFI, "[", 2)
    names(iterations) <- time_span
    iterations <- do.call("rbind", iterations)

    convergence <- sapply(RKFI, "[", 3)
    names(convergence) <- time_span
    convergence <- do.call("rbind", convergence)

    diversification <- sapply(RKFI, "[", 4)
    diversification <- do.call("rbind.data.frame", diversification)

    ubiquity <- sapply(RKFI, "[", 5)
    ubiquity <- do.call("rbind.data.frame", ubiquity)

    RKFI <- sapply(RKFI, "[", 1)
    RKFI <- do.call("rbind.data.frame", RKFI)

    class(RKFI) <- c('data.frame', 'reks_fitness_tccgp')
    attr(RKFI, 'diversification') <- diversification
    attr(RKFI, 'ubiquity') <- ubiquity
    attr(RKFI, "iterations") <- iterations
    attr(RKFI, "convergence") <- convergence
    if (binary_mode == 'RTA') {
        attr(RKFI, "binary_mode") <- 'RTA'
    }
    if (binary_mode == 'RCA') {
        attr(RKFI, "binary_mode") <- 'RCA'
    }
    if (binary_mode == 'simple') {
        attr(RKFI, "binary_mode") <- 'simple'
    }
    if (binary_mode == "higher_quartiles") {
        attr(RKFI, "binary_mode") <- 'higher_quartiles'
    }
    if (binary_mode == "higher_quartiles_kng") {
        attr(RKFI, "binary_mode") <- 'higher_quartiles_kng'
    }
    attr(RKFI, 'geo_dim') <- geo_dim
    attr(RKFI, 'kng_dim') <- kng_dim
    attr(RKFI, 'time_dim') <- time_dim
    attr(RKFI, 'measure') <- measure

    # TODO
    # It seems there's some problem about the use of memory, but I don't
    #  know if it's the right way to solve the problem
    gc()

    return(RKFI)
}
