#' @name
#' coherence_panel
#'
#' @title
#' Regional Coherence Index for panel data
#'
#' @description
#' The function computes the so called Coherence index.
#'
#' @details
#' The function computes the so called Coherence index.
#' It assumes that the "universe" from which you derive the distribution of
#' reference is composed by all the information provided in the database, and
#' only this.
#' You can use it both on a panel data set (if you identify also a column with
#' a temporal indication of the observations) or on cross-section data (by
#' leaveing the parameter just said as NULL).
#' See:
#' \itemize{
#' \item{Teece, Rumelt, Dosi and Winter (1994) 'Understanding Corporate
#' Coherence: Theory and Evidence', \emph{Journal of Economic Behavior \&
#' Organization}, 23, 1-30;}
#' \item{Nesta and Saviotti (2005) 'Coherence of the Knowledge Base and the
#' Firm's Innovative Performance: Evidence from the U.S. Pharmaceutical
#' Industry', \emph{Journal of Industrial Economics}, 53, 123-142;}
#' \item{Nesta and Saviotti (2006) 'Firm Knowledge and Market Value in
#' Biotechnology', \emph{Industrial and Corporate Change}, 15, 625-652;}
#' \item{Quatraro (2010) 'Knowledge Coherence, Variety and Economic Growth:
#' Manufacturing Evidence from Italian Regions', \emph{Research Policy}, 39,
#' 1289-1302.}
#'
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
#' @return A data.frame with the Coherence Index of each geographical area
#' in each time step of analysis.
#' @examples
#' RCI <- coherence(data = df, time_dim = year, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

coherence_panel <- function(data, geo_dim, kng_dim, kng_nbr, time_dim) {

    # Preliminary operations, checks and transformations ---------------

    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop(paste0("Package \"reshape2\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
        stop(paste0("Package \"plyr\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    data <- as.data.frame(data)

    time_dim <- deparse(substitute(time_dim))
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))

    data <- data[complete.cases(data[, c(time_dim, geo_dim,
                                         kng_dim, kng_nbr)]), ]

    # Null model (Hypergeometric distribution) ---------------

    C <- .get_biadj_matrix(data, geo_dim, kng_dim)

    J <- t(C) %*% C

    n <- colSums(C)

    mu <- n %*% t(n)
    mu <- mu/nrow(C)
    rownames(mu) <- colnames(mu)

    s2 <- (1 - n/nrow(C)) %*% t((nrow(C) - n)/(nrow(C) - 1))
    s2 <- mu * s2
    rownames(s2) <- colnames(s2)

    t <- (J - mu)/sqrt(s2)

    diag(t) <- 0

    null_model <- list()
    null_model$t <- t
    class(null_model) <- "reks_null_model"

    # Coherence index for panel data ---------------

    RCI <- lapply(unique(data[, time_dim]), function(t) {
        data_subset <- data[which(data[, time_dim] == t), ]
        R <- coherence(data_subset, geo_dim, kng_dim, kng_nbr,
                       null_model = null_model, names_as_strings = TRUE)
        colnames(R) <- c(geo_dim, as.character(t))
        return(R)
    })

    RCI <- plyr::join_all(RCI, geo_dim)
    RCI <- reshape2::melt(RCI)
    colnames(RCI) <- c(geo_dim, time_dim, "Coherence")

    return(RCI)
}
