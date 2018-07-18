#' @name
#' fitness_tccgp
#'
#' @title
#' Regional Knowledge Fitness Index a là Tacchella, Cristelli,
#' Caldarelli, Gabrielli and Pietronero
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
#' @param names_as_strings If TRUE the function assumes you have written the
#' name of the columns above said as strings, otherwise they are assumed to
#' be symbols. Default is FALSE.
#' @return A data.frame with the Fitness Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RKFI <- fitness(data = df, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

fitness_tccgp <- function(data, geo_dim, kng_dim, kng_nbr,
                          binary_mode = "RTA", names_as_strings = FALSE) {
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

    # This is not needed for the algorithm, but still it can be useful to have
    #  this information stored for future purpouses.
    du <- .get_du(mm)

    RKFI <- rep(1, nrow(mm))
    RKCI <- rep(1, ncol(mm))
    i <- 0
    while (TRUE) {
        RKFI1 <- rowSums(t(t(mm) * RKCI))
        RKCI1 <- 1 / rowSums(t(mm) / RKFI)

        # Normalisation needed to avoid possible divergences
        #  due to the hyperbolic nature of the second equation
        RKFI1 <- RKFI1 / mean(RKFI1)
        RKCI1 <- RKCI1 / mean(RKCI1)

        # TODO
        # This is an arbitrary choice of mine and it is not in the
        #  original papers. Maybe it could be better and closer to the
        #  original sources to use 20 recursions. Another thing to check
        #  and decide is what happens in the function if the algorithm
        #  does not converge
        if (all((RKFI - RKFI1) < 0.0000000001) &
            all((RKCI - RKCI1) < 0.0000000001)) {
            RKFI <- RKFI1
            RKCI <- RKCI1
            convergence <- TRUE
            break()
        }
        if (i >= 200) {
            RKFI <- rep(as.numeric(NA), nrow(mm))
            RKCI <- rep(as.numeric(NA), ncol(mm))
            names(RKFI) <- names(RKFI1)
            names(RKCI) <- names(RKCI1)
            convergence <- FALSE
            warning(paste0('The algorithm failed to converge.\n',
                           'Maybe your matrix is not triangular ',
                           'as expected.\nYou can check it using ',
                           'plot_biadj_matrix()'))
            break()
        }

        RKFI <- RKFI1
        RKCI <- RKCI1

        i <- i + 1
    }

    RKFI <- cbind.data.frame(names(RKFI), RKFI)
    measure <- "Fitness"
    colnames(RKFI) <- c(geo_dim, measure)

    class(RKFI) <- c('data.frame', 'reks_fitness_tccgp')
    attr(RKFI, 'diversification') <- du$diversification
    attr(RKFI, 'ubiquity') <- du$ubiquity
    attr(RKFI, "iterations") <- i
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
    attr(RKFI, 'measure') <- measure

    # TODO
    # It seems there's some problem about the use of memory, but I don't
    #  know if it's the right way to solve the problem
    gc()

    return(RKFI)
}
