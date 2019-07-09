#' @name
#' complexity
#'
#' @title
#' Regional (Technological) Complexity Index
#'
#' @description
#' The function computes the complexity index of each region, in line with the
#' methodology proposed by Hidalgo, Hausmann and coauthors.
#'
#' @references
#' Hidalgo and Hausmann (2009) ``The Building Blocks of Economic
#' Complexity'', \emph{PNAS}, 106, 10570--10575;
#'
#' Hausmann, Hidalgo, Bustos, Coscia, Chung, Jimenez, Simoes, and Yildirim
#' (2014) \emph{The Atlas of Economic Complexity}, The MIT Press, 1st ed. 2011.
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
#' @importFrom stats cor
#'
#' @param occt Contingency table (i.e., occurrence table or incidence
#' matrix) on which you want to compute the indices. It can be a 2D array,
#' in which the first dimension represents the units of analysis (like firms,
#' regions, or countries), and the second dimension represents the
#' events or characteristics of interest (like the classes of the patents
#' produced by the regions, or the sectors in which the workers belongs).
#' Lastly, the values in each cell represents the occurrences of each unit-event
#' pair. Moreover, you can use also a 3D array if you like, in which the third
#' dimension represents the time. The object is expected to be of "table" class.
#' @param rta If TRUE (default) it uses the Revealed Technological Advantages
#' (RTA) of the original data
#' @param binary If TRUE (default) it dichotomize the RTA matrix
#' (can be used only together with rta=TRUE).
#' @param scale If TRUE (default) the output is standardised (mean = 0;
#' sd = 1).
#' @param which It can be one of "rci" (default), "tci", "both". The first
#' returns the complexity of each region; the second returns the complexity of
#' each technological domain; and the third returns both the indices.
#'
#' @return A data.frame with the Complexity Index of each region and/or of each
#' technological domain. If a 3D array is provided as input, it returns the
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
#' CX <- complexity(octab)
#'
#' @export

complexity <- function(occt,
                       rta = TRUE,
                       binary = TRUE,
                       scale = TRUE,
                       which = "rci") {

    #-- Preliminary steps and checks
    info <- data_info(occt)

    #-- Main Function
    if (info$n_dims == 3) {
        index <- apply(occt, 3, complexity,
                    rta, binary, scale, which)
        if (which == "both") {
            RCI <- lapply(index, function(x) {
                unique(x[, c(1, 3)])
            })
            TCI <- lapply(index, function(x) {
                unique(x[, c(2, 4)])
            })
            index <- list(RCI = RCI,
                          TCI = TCI)
        }
    } else {
        occt <- remove_zeros(occt)
        if (isTRUE(rta)) {
            occt <- rta(occt, binary)
        }
        diversity <- Matrix::rowSums(occt)
        ubiquity <- Matrix::colSums(occt)
        mcp1 <- occt / diversity
        mcp2 <- Matrix::t(Matrix::t(occt) / ubiquity)
        # Mcc <- Matrix::tcrossprod(mcp1, mcp2)
        Mpp <- Matrix::crossprod(mcp2, mcp1)
        if (!all(round(Matrix::rowSums(Mpp)) == 1)) {
            stop(paste("The matrix is not row-stochastic.\n",
                       "It is not possible to compute the measure."))
        }
        if (round(Re(as.complex(eigen(Mpp)$value[1]))) != 1) {
            stop(paste("The first eigen-value is different from 1.",
                       "It is not possible to compute the measure."))
        }
        TCI <- eigen(Mpp)$vectors
        if (dim(TCI)[2] >= 2) {
            TCI <- TCI[, 2]
            TCI <- Re(as.complex(TCI))
            RCI <- as.vector(mcp1 %*% TCI)
            if (cor(RCI, diversity,
                    use = "pairwise.complete.obs",
                    method = "spearman") < 0) {
                RCI <- -RCI
                TCI <- -TCI
            }
        } else {
            RCI <- NA
            TCI <- NA
        }
        if (scale == TRUE) {
            RCI <- scale(RCI)
            # Note: PCI values are normalized using ECI mean and st dev
            #  to preserve property that ECI = (mean PCI of products
            #  for which Mcp = 1)
            TCI <- (TCI - attr(RCI, "scaled:center"))
            TCI <- TCI / attr(RCI, "scaled:scale")
        }
        names(RCI) <- rownames(occt)
        names(TCI) <- colnames(occt)
        index <- switch(which,
            rci = RCI,
            tci = TCI,
            both = list(RCI = RCI,
                        TCI = TCI))
    }

    #-- Final steps
    index <- switch(which,
                    rci = wide_to_long(index,
                                       info$n_dims,
                                       info$dim_nms[info$nd[[1]]],
                                       "RCI"),
                    tci = wide_to_long(index,
                                       info$n_dims,
                                       info$dim_nms[info$nd[[2]]],
                                       "TCI"),
                    both = merge(wide_to_long(index$RCI,
                                              info$n_dims,
                                              info$dim_nms[info$nd[[1]]],
                                              "RCI"),
                                 wide_to_long(index$TCI,
                                              info$n_dims,
                                              info$dim_nms[info$nd[[2]]],
                                              "TCI"),
                                 all = TRUE)[, c(info$n_dims - 1,
                                                 info$n_dims + 1,
                                                 info$n_dims - 2,
                                                 info$n_dims,
                                                 info$n_dims + 2)]
    )
    return(index)
}
