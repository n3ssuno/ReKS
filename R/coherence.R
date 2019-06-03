#' @name
#' coherence
#'
#' @title
#' Coherence
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
#'
#' @references
#' Teece, Rumelt, Dosi and Winter (1994) ``Understanding Corporate
#' Coherence: Theory and Evidence'', \emph{Journal of Economic Behavior &
#' Organization}, 23, 1--30;
#'
#' Nesta and Saviotti (2005) ``Coherence of the Knowledge Base and the
#' Firm's Innovative Performance: Evidence from the {U.S.}~Pharmaceutical
#' Industry'', \emph{Journal of Industrial Economics}, 53, 123--142;
#'
#' Nesta and Saviotti (2006) ``Firm Knowledge and Market Value in
#' Biotechnology'', \emph{Industrial and Corporate Change}, 15, 625--652;
#'
#' Nesta (2008) ``Knowledge and Productivity in the World's Largest
#' Manufacturing Corporations'', \emph{Journal of Economic Behavior \&
#' Organization}, 67, 886--902,
#'
#' Bottazzi and Pirino (2010) ``Measuring Industry Relatedness and
#' Corporate Coherence'', \emph{SSRN Electronic Journal}, 11, 1--24;
#'
#' Quatraro (2010) ``Knowledge Coherence, Variety and Economic Growth:
#' Manufacturing Evidence from Italian Regions'', \emph{Research Policy}, 39,
#' 1289-1302.
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
#' @param relatedness_mtx Matrix of similarity between the technological
#' classes. The ReKS package provides some function to buid it: see
#' \code{\link{relatedness}} and \code{\link{proximity}}.
#'
#' @return A data.frame with the Coherence Index of each geographical area.
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
#' rel_m <- relatedness(octab, output_statistic = "p")
#' RCH <- coherence(octab, rel_m)
#'
#' @export

coherence <- function(occt, relatedness_mtx) {

    #-- Preliminary steps and checks
    info <- data_info(occt)

    #-- Main function
    if (info$n_dims == 3) {
        C <- apply(as.table(occt), 3, coherence, relatedness_mtx)
    } else {
        # Preliminary operations, checks and transformations
        occt <- Matrix::Matrix(occt)
        oc_mtx_names <- colnames(occt)
        rl_mtx_names <- colnames(relatedness_mtx)
        if (dim(occt)[[2]] != sum(dim(relatedness_mtx)) / 2) {
            names_tbr <- setdiff(rl_mtx_names, oc_mtx_names)
            if (length(names_tbr) != 0)
                relatedness_mtx <- relatedness_mtx[
                    -which(rownames(relatedness_mtx) %in% names_tbr),
                    -which(colnames(relatedness_mtx) %in% names_tbr)]
        }
        if (dim(occt)[[2]] != sum(dim(relatedness_mtx)) / 2)
            stop(paste('There is some problem, because the two matrices',
                       'considered have a different number of',
                       'columns.'), call. = FALSE)
        rl_mtx_names <- colnames(relatedness_mtx)
        if (any(oc_mtx_names != rl_mtx_names)) {
            oc_mtx_names <- oc_mtx_names[, order(colnames(oc_mtx_names))]
            relatedness_mtx <- relatedness_mtx[order(rownames(relatedness_mtx)),
                                               order(colnames(relatedness_mtx))]
        }
        if (any(oc_mtx_names != rl_mtx_names))
            stop(paste('There is some problem, because the there is no perfect',
                       'correspondence between the column names of the two',
                       'matrices considered.'), call. = FALSE)
        # Waighted Average Relatedness
        ones <- !Matrix::diag(TRUE,
                              nrow = nrow(relatedness_mtx),
                              ncol = nrow(relatedness_mtx))
        mode(ones) <- "integer"
        WAR_num <- Matrix::tcrossprod(occt, relatedness_mtx)
        WAR_den <- Matrix::tcrossprod(occt, ones)
        WAR <- WAR_num / WAR_den
        # Coherence
        C_num <- WAR * occt
        C_den <- Matrix::rowSums(occt)
        C <- Matrix::rowSums(C_num / C_den)
        C[which(is.nan(C))] <- 0
    }

    #-- Final steps
    C <- wide_to_long(C, info$n_dims,
                      info$dim_nms[info$nd[[1]]], "Coherence")
    return(C)
}

