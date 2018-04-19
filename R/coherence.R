#' @name
#' coherence
#'
#' @title
#' Regional Coherence Index
#' @description
#' a
#' @details
#' See:
#' \itemize{
#' \item{Teece, Rumelt, Dosi and Winter (1994) 'Understanding Corporate
#'        Coherence: Theory and Evidence', \emph{Journal of Economic
#'        Behavior \& Organization}, 23, 1-30;}
#' \item{Nesta and Saviotti (2005) 'Coherence of the Knowledge Base and the
#'        Firm's Innovative Performance: Evidence from the U.S. Pharmaceutical
#'        Industry', \emph{Journal of Industrial Economics}, 53, 123-142;}
#' \item{Nesta and Saviotti (2006) 'Firm Knowledge and Market Value in
#'        Biotechnology', \emph{Industrial and Corporate Change}, 15, 625-652;}
#' \item{Quatraro (2010) 'Knowledge Coherence, Variety and Economic Growth:
#'        Manufacturing Evidence from Italian Regions', \emph{Research Policy},
#'        39, 1289-1302.}
#' }
#' @param data It is expected to be a data.frame in 'long' format.
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
#' It is an optional parameter
#' @return A data.frame with the Coherence Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RCI <- Coherence(data = regpat, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

Coherence <- function(data, geo_dim, kng_dim, kng_nbr, time_dim = NULL) {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop(paste0("Package \"reshape2\" needed for this function to work. ",
            "Please install it."), call. = FALSE)
    }
    library(reshape2)
    if (!requireNamespace("plyr", quietly = TRUE)) {
        stop(paste0("Package \"plyr\" needed for this function to work. ",
            "Please install it."), call. = FALSE)
    }
    library(plyr)

    data <- as.data.frame(data)

    if (is.null(time_dim)) {
        data$AddedTimeDim <- 1
        time_dim <- "AddedTimeDim"
        time_dim_added <- TRUE
    }

    if (!is.character(time_dim)) {
        time_dim <- deparse(substitute(time_dim))
    }
    if (!is.character(geo_dim)) {
        geo_dim <- deparse(substitute(geo_dim))
    }
    if (!is.character(kng_dim)) {
        kng_dim <- deparse(substitute(kng_dim))
    }
    if (!is.character(kng_nbr)) {
        kng_nbr <- deparse(substitute(kng_nbr))
    }

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    C <- unique(data[, c(geo_dim, kng_dim)])
    C$C <- 1
    C <- dcast(C, formula(paste(geo_dim, "~", kng_dim)), value.var = "C")
    C[is.na(C)] <- 0
    gnames <- C[, 1]
    C <- as.matrix(C[, -1])
    rownames(C) <- gnames

    J <- t(C) %*% C

    n <- colSums(C)

    mu <- n %*% t(n)
    mu <- mu/nrow(C)
    rownames(mu) <- colnames(mu)

    s2 <- (1 - n/nrow(C)) %*% t((nrow(C) - n)/(nrow(C) - 1))
    s2 <- mu * s2
    rownames(s2) <- colnames(s2)

    t <- (J - mu)/sqrt(s2)

    t1 <- t
    diag(t1) <- 0

    ones <- matrix(1, nrow = nrow(t), ncol = nrow(t))
    ones1 <- ones
    diag(ones1) <- 0

    RCI <- lapply(unique(data[, time_dim]), function(t) {
        ee <- dcast(data[which(data[, time_dim] == t), ],
                    formula(paste(geo_dim, "~", kng_dim)),
                    value.var = kng_nbr, fun.aggregate = sum)
        gnames <- ee[, 1]
        ee <- as.matrix(ee[, -1])
        rownames(ee) <- gnames
        en <- colnames(ee)
        nc <- setdiff(colnames(t1), colnames(ee))
        en <- c(en, nc)
        for (i in nc) {
            ee <- cbind(ee, 0)
        }
        colnames(ee) <- en
        col.order <- colnames(t1)
        ee <- ee[, col.order]

        WAR <- ee %*% t1/ee %*% ones1
        R <- rowSums(WAR * ee/rowSums(ee))
        R[is.nan(R)] <- 0
        R <- as.data.frame(R)
        R <- cbind(rownames(R), R)
        colnames(R) <- c(geo_dim, as.character(t))
        return(R)
    })
    names(RCI) <- unique(data[, time_dim])
    RCI <- join_all(RCI, geo_dim)
    RCI <- melt(RCI)
    colnames(RCI) <- c(geo_dim, time_dim, "Coherence.Index")
    if (time_dim_added) {
        RCI <- RCI[, c(geo_dim, "Coherence.Index")]
    }
    return(RCI)
}
