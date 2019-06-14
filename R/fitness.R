#' @name
#' fitness
#'
#' @title
#' Regional Fitness and/or Technological Complexity
#'
#' @description
#' The function computes the fitness index of each region (i.e.
#' its \emph{competitiveness}), in line with the methodology proposed by
#' Pietronero and coauthors.
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
#' @param occt Contingency table (i.e., occurrence table or incidence
#' matrix) on which you want to compute the indices. It can be a 2D array,
#' in which the first dimension represents the units of analysis (like firms,
#' regions, or countries), and the second dimension represents the
#' events or characteristics of interest (like the classes of the patens
#' produced by the regions, or the sectors in which the workers belongs).
#' Lastly, the values in each cell represents the occurrences of each unit-event
#' pair. Moreover, you can use also a 3D array if you like, in which the third
#' dimension represents the time. The object is expected to be of "table" class.
#' @param rta If TRUE (default) it uses the Revealed Technological Advantages
#' (RTA) of the original data
#' @param binary If TRUE (default) it binarises the RTA matrix
#' (can be used only together with rta=TRUE).
#' @param which It can be one of "fitness" (default), "complexity", "both".
#' The first returns the fitness of each region; the second returns the
#' complexity of each technological domain; and the third returns both the
#' indices.
#'
#' @return A data.frame with the Fitness of each region and/or the Complexity of
#' each technological domain. If a 3D array is provided as input, it returns the
#' full panel data.frame.
#'
#' @examples
#' geo <- paste0("R", 10:50)
#' tech <- paste0("T", 10:90)
#' time <- 1:5
#' dat <- expand.grid(geo, tech, time)
#' colnames(dat) <- c("geo", "tech", "time")
#' set.seed(1)
#' dat$nPat <- sample(1:200, nrow(dat), replace = TRUE)
#' octab <- xtabs(nPat ~ geo + tech + time, dat)
#' octab[sample(1:length(octab), length(octab)/2)] <- 0
#' FX <- fitness(octab)
#' attr(FX, "convergence")
#'
#' @export

# @importFrom methods as is

fitness <- function(occt,
                    rta = TRUE,
                    binary = TRUE,
                    which = "fitness") {

    #-- Preliminary steps and checks
    info <- data_info(occt)

    #-- Main function
    if (info$n_dims == 3) {
        index <- apply(occt, 3, fitness, rta, binary, which)
        i <- sapply(index, attr, "iterations")
        c <- sapply(index, attr, "convergence")
        if (which == "both") {
            Fx <- lapply(index, function(x) {
                unique(x[, c(1, 3)])
            })
            Cx <- lapply(index, function(x) {
                unique(x[, c(2, 4)])
            })
            index <- list(Fx = Fx,
                          Cx = Cx)
        }
    } else {
        occt <- remove_zeros(occt)
        spm <- any(methods::is(occt) == "sparseMatrix")
        if (isTRUE(rta)) {
            occt <- rta(occt, binary)
        }
        Fx <- rep(1, nrow(occt))
        Cx <- rep(1, ncol(occt))
        if (spm) {
            Fx <- Matrix::c.sparseVector(Fx)
            Cx <- Matrix::c.sparseVector(Cx)
        }
        i <- 0
        while (TRUE) {
            Fx1 <- Matrix::t(occt) * Cx
            Fx1 <- Matrix::rowSums(Matrix::t(Fx1))
            Cx1 <- 1 / Matrix::rowSums(Matrix::t(occt / Fx))
            # Normalisation needed to avoid possible divergences
            #  due to the hyperbolic nature of the second equation
            Fx1 <- Fx1 / mean(Fx1)
            Cx1 <- Cx1 / mean(Cx1)
            # TODO
            # This is an arbitrary choice of mine and it is not in the
            #  original papers. Maybe it could be better and closer to the
            #  original sources to use 20 recursions. Another thing to check
            #  and decide is what happens in the function if the algorithm
            #  does not converge
            if (all(Fx - Fx1 < 0.0000000001) &
                all(Cx - Cx1 < 0.0000000001)) {
                Fx <- Fx1
                Cx <- Cx1
                c <- TRUE
                break
            }
            if (i >= 200) {
                Fx <- rep(as.numeric(NA), nrow(occt))
                Cx <- rep(as.numeric(NA), ncol(occt))
                names(Fx) <- names(Fx1)
                names(Cx) <- names(Cx1)
                c <- FALSE
                break
            }
            Fx <- Fx1
            Cx <- Cx1
            i <- i + 1
        }
        index <- switch(which,
                        fitness = Fx,
                        complexity = Cx,
                        both = list(Fx = Fx,
                                    Cx = Cx))
    }

    #-- Final steps
    index <- switch(which,
                    fitness = wide_to_long(index,
                                           info$n_dims,
                                           info$dim_nms[info$nd[[1]]],
                                           "fitness"),
                    complexity = wide_to_long(index,
                                              info$n_dims,
                                              info$dim_nms[info$nd[[2]]],
                                              "complexity"),
                    both = merge(wide_to_long(index$Fx,
                                              info$n_dims,
                                              info$dim_nms[info$nd[[1]]],
                                              "fitness"),
                                 wide_to_long(index$Cx,
                                              info$n_dims,
                                              info$dim_nms[info$nd[[2]]],
                                              "complexity"),
                                 all = TRUE)[, c(info$n_dims - 1,
                                                 info$n_dims + 1,
                                                 info$n_dims - 2,
                                                 info$n_dims,
                                                 info$n_dims + 2)]
    )
    attr(index, "iterations") <- i
    attr(index, "convergence") <- c
    return(index)
}
