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

fitness_tccgp <- function(occurrence_mtx,
                          rta = TRUE, binary = TRUE, scale = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop(paste0('Package \"Matrix\" needed for this function to work. ',
                    'Please install it.'), call. = FALSE)

    geo_dim <- attr(occurrence_mtx, 'geo_dim')
    kng_dim <- attr(occurrence_mtx, 'kng_dim')
    if (is.list(occurrence_mtx))
        time_dim <- attr(occurrence_mtx, 'time_dim')
    measure <- "Fitness"

    fitness_tccgp_crossSection <- function(occurrence_mtx) {
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

        if (isTRUE(rta))
            occurrence_mtx <- rta(occurrence_mtx, binary = binary)
        if (isTRUE(binary))
            occurrence_mtx <- Matrix::Matrix(ifelse(occurrence_mtx > 0, 1, 0),
                                             nrow = nrow(occurrence_mtx))

        # This is not needed for the algorithm, but still it can be useful to have
        #  this information stored for future purpouses.
        du <- ReKS:::.get_du(Matrix::as.matrix(occurrence_mtx))
        #TODO

        RKFI <- as(rep(1, nrow(occurrence_mtx)), "sparseVector")
        RKCI <- as(rep(1, ncol(occurrence_mtx)), "sparseVector")
        i <- 0
        while (TRUE) {
            RKFI1 <- Matrix::t(occurrence_mtx) * RKCI
            RKFI1 <- Matrix::rowSums(Matrix::t(RKFI1))
            RKCI1 <- 1 / Matrix::rowSums(Matrix::t(occurrence_mtx / RKFI))

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
                RKFI <- rep(as.numeric(NA), nrow(occurrence_mtx))
                RKCI <- rep(as.numeric(NA), ncol(occurrence_mtx))
                names(RKFI) <- names(RKFI1)
                names(RKCI) <- names(RKCI1)
                convergence <- FALSE
                warning(paste0('The algorithm failed to converge.\n',
                               'Maybe your matrix is not triangular ',
                               'as expected.\nYou can check it using ',
                               'image(occurrence_mtx, useAbs = FALSE)'))
                break()
            }

            RKFI <- RKFI1
            RKCI <- RKCI1

            i <- i + 1
        }

        RKFI <- cbind.data.frame(rnms,
                                 RKFI)
        colnames(RKFI) <- c(geo_dim, measure)

        if (scale == TRUE) {
            RKFI[, measure] <- scale(as.numeric(RKFI[, measure]))
            # warning('The values of the index have been standardised.')
        }

        # TODO
        # It seems there's some problem about the use of memory, but I don't
        #  know if it's the right way to solve the problem
        gc()

        attr(RKFI, "iterations") <- i
        attr(RKFI, "convergence") <- convergence

        attr(RKFI, 'diversification') <- du$diversification
        attr(RKFI, 'ubiquity') <- du$ubiquity

        return(RKFI)
    }

    fitness_tccgp_panel <- function(occurrence_mtx) {
        time_span <- names(occurrence_mtx)
        RKFI <- lapply(as.character(time_span),
                       function(y) {
                           Fx <- fitness_tccgp_crossSection(occurrence_mtx[[y]])

                           iterations <- attr(Fx, "iterations")
                           convergence <- attr(Fx, "convergence")

                           # diversification <- attr(Fx, 'diversification')
                           # diversification <-
                           #     cbind.data.frame(geo_dim = names(diversification),
                           #                      time_dim = y,
                           #                      diversification)
                           # ubiquity <- attr(Fx, 'ubiquity')
                           # ubiquity <-
                           #     cbind.data.frame(kng_dim = names(ubiquity),
                           #                      time_dim = y,
                           #                      ubiquity)

                           return(list(Fx,
                                       iterations, convergence))
                                       # diversification, ubiquity))
                       })

        iterations <- sapply(RKFI, "[", 2)
        names(iterations) <- time_span
        iterations <- do.call("rbind", iterations)

        convergence <- sapply(RKFI, "[", 3)
        names(convergence) <- time_span
        convergence <- do.call("rbind", convergence)

        # diversification <- sapply(RKFI, "[", 4)
        # diversification <- do.call("rbind.data.frame", diversification)
        #
        # ubiquity <- sapply(RKFI, "[", 5)
        # ubiquity <- do.call("rbind.data.frame", ubiquity)

        RKFI <- sapply(RKFI, "[", 1)
        yrs <- unlist(mapply(rep, time_span, lapply(RKFI, nrow)))
        RKFI <- do.call("rbind.data.frame", RKFI)
        RKFI <- cbind.data.frame(RKFI, yrs)
        RKFI <- RKFI[, c(1, 3, 2)]
        colnames(RKFI) <- c(geo_dim, time_dim, measure)

        # attr(RKFI, 'diversification') <- diversification
        # attr(RKFI, 'ubiquity') <- ubiquity
        attr(RKFI, "iterations") <- iterations
        attr(RKFI, "convergence") <- convergence

        return(RKFI)
    }

    fntn <- ifelse(is.list(occurrence_mtx),
                   "fitness_tccgp_panel",
                   "fitness_tccgp_crossSection")
    Fx <- do.call(fntn, list(occurrence_mtx))

    # final steps

    # class(R) <- c("reks_fitness_tccgp", "data.frame")
    attr(Fx, 'geo_dim') <- geo_dim
    attr(Fx, 'kng_dim') <- kng_dim
    if (is.list(occurrence_mtx))
        attr(Fx, 'time_dim') <- time_dim
    attr(Fx, 'measure') <- measure
    # attr(Fx, 'diversity') <- du$diversification
    # attr(Fx, 'ubiquity') <- du$ubiquity
    attr(Fx, 'standardised') <- scale
    attr(Fx, "RTA") <- rta
    attr(Fx, "binary") <- binary
    # attr(Fx, "iterations") <- i
    # attr(Fx, "convergence") <- convergence

    return(Fx)
}
