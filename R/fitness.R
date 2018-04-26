#' @name
#' fitness
#'
#' @title
#' Regional Knowledge Fitness Index
#'
#' @description
#'
#' @details
#'
#' See:
#' \itemize{
#' \item{Tacchella, Cristelli, Caldarelli, Gabrielli and Pietronero (2012)
#' "A New Metrics for Countries' Fitness and Products' Complexity",
#' \emph{Scientific Reports}, 2, 1--7;}
#' \item{Tacchella, Cristelli, Caldarelli, Gabrielli and Pietronero (2013)
#' "Economic Complexity: Conceptual Grounding of a New Metrics for Global
#' Competitiveness", \emph{Journal of Economic Dynamics and Control}, 37,
#' 1683--1691;}
#' \item{Cristelli, Gabrielli, Tacchella, Caldarelli and Pietronero (2013)
#' "Measuring the Intangibles: A Metrics for the Economic Complexity of
#' Countries and Products", \emph{PLOS ONE}, 8, e70726;}
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
#' assuming that some are more "ubiquitous" than others). Warning: look
#' carefully at the plot_biadj_matrix() output because I personally observed
#' that with RTA applied to patent data you do not have the triangular structure
#' observed in the trade data.
#' @return A data.frame with the Fitness Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RKFI <- fitness(data = df, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

fitness <- function(data, geo_dim, kng_dim, kng_nbr, time_dim = NULL,
                    binary_mode = "RTA") {
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

    RKFI <- lapply(unique(data[, time_dim]), function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]
        mm <- .get_biadj_matrix(data_subset, geo_dim, kng_dim, kng_nbr,
                                binary_mode)

        if (any(rowSums(mm) == 0)) {
            mm <- mm[-which(rowSums(mm) == 0), ]
        }
        if (any(colSums(mm) == 0)) {
            mm <- mm[, -which(colSums(mm) == 0)]
        }

        ff <- rep(1, nrow(mm))
        cc <- rep(1, ncol(mm))
        i <- 0
        while (TRUE) {
            ff1 <- rowSums(t(t(mm) * cc))
            cc1 <- 1 / rowSums(t(mm) / ff)

            # Normalisation needed to avoid possible divergences
            #  due to the hyperbolic nature of the second equation
            ff1 <- ff1 / mean(ff1)
            cc1 <- cc1 / mean(cc1)

            # TODO
            # This is an arbitrary choice of mine and it is not in the
            #  original papers. Maybe it could be better and closer to the
            #  original sources to use 20 recursions. Another thing to check
            #  and decide is what happens in the function if the algorithm
            #  does not converge
            if (all((ff - ff1) < 0.0000000001) &
                all((cc - cc1) < 0.0000000001)) {
                ff <- ff1
                cc <- cc1
                # print(paste0(t, ": ", i))
                break()
            }

            ff <- ff1
            cc <- cc1

            i <- i + 1
        }

        ff <- cbind.data.frame(names(ff), ff)
        colnames(ff) <- c(geo_dim, as.character(t))

        return(list(ff, i))
        # return(ff)
    })

    iterations <- sapply(RKFI, "[", 2)
    iterations <- unlist(iterations)
    if (!time_dim_added) {
        names(iterations) <- unique(data[, time_dim])
    }

    RKFI <- sapply(RKFI, "[", 1)
    names(RKFI) <- unique(data[, time_dim])
    RKFI <- plyr::join_all(RKFI, geo_dim)
    RKFI <- reshape2::melt(RKFI)
    colnames(RKFI) <- c(geo_dim, time_dim, "Fitness")
    if (time_dim_added) {
        RKFI <- RKFI[, c(geo_dim, "Fitness")]
    }

    class(RKFI) <- c('data.frame', 'rks_fitness')
    # TODO
    # In that way it returns only the last year
    # attr(RKFI, 'diversification') <- du$diversification
    # attr(RKFI, 'ubiquity') <- du$ubiquity
    attr(RKFI, "iterations") <- iterations
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

    # TODO
    # It seems there's some problem about the use of memory, but I don't
    #  know if it's the right way to solve the problem
    gc()


    return(RKFI)
}
