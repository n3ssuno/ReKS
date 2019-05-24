#' @name
#' shannonEntropy
#'
#' @title
#' Shannon Entropy and Theil entropy decomposition
#'
#' @description
#' The function returns the Shannon information entropy of a given table.
#' Optionally, it can return also the values for the entropy decomposed in
#' the two constituents of within- and between-groups entropy.
#'
#' @details
#' The function returns the Shannon information entropy of a given table.
#' It computes internally the frequency of each observation and, if not
#' specified differently, uses "bits" as units of measure.
#' The function can compute either the entropy of a cross-section, or of
#' a panel.
#' If a vector of groups is provided, the function returns also the values
#' for the entropy decomposed in the two constituents of within- and
#' between-groups entropy, following the Theil's decomposition.
#'
#' @references
#' Theil (1972) \emph{Statistical Decomposition Analysis}, North-Holland.
#'
#' Zadjenweber (1972) ``Une Application de la Théorie de l'Information
#' à l'Économie: La Mesure de la Concentration'', \emph{Revue d'Economie
#' Politique}, 82, 486--510.
#'
#' Frenken (2007) ``Entropy Statistics and Information Theory'', in
#' Hanusch and Pyka (Eds.) \emph{Elgar Companion to Neo-Schumpeterian
#' Economics}, Edward Elgar.
#'
#' Stirling (2007) ``A General Framework for Analysing Diversity in Science,
#' Technology and Society'', \emph{Journal of The Royal Society Interface}, 4,
#' 707--719.
#'
#' Quatraro (2010) ``Knowledge Coherence, Variety and Economic Growth:
#' Manufacturing Evidence from Italian Regions'', \emph{Research Policy}, 39,
#' 1289-1302;
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
#' @param log_base Base of the logarithm used to compute the entropy (i.e.,
#' its unit of measure). The default value is 2. You can use any number you
#' like (to use the natural logaritihm, write exp(1)).
#' @param Theil_decomp A list of groups . It should be of the same length
#' of the second dimension of the occt provided.
#'
#' @return data.table with the (total) entropy and (eventually) its
#' decomposition in between- and within-group components.
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
#' grps <- substr(attr(octab, "dimnames")$tech, 2, 2)
#' ETP <- shannonEntropy(octab,
#'                       log_base = exp(1),
#'                       Theil_decomp = grps)

shannonEntropy <- function(occt,
                           log_base = 2,
                           Theil_decomp = NULL) {

    #-- Preliminary steps and checks
    info <- data_info(occt)

    #-- Internal Functions
    entropy <- function(d, log_base) {
        d <- d[which(d > 0)]
        f <- d / sum(d)
        e <- -sum(f * log(f, log_base))
        e[!is.finite(e)] <- 0
        return(e)
    }

    #-- Main Function
    if (info$n_dims == 3) {
        etp <- apply(occt, 3, shannonEntropy, log_base, Theil_decomp)
        if(!is.null(Theil_decomp)) {
            betp <- lapply(etp, `[`, c(info$dim_nms[1], "entropy_between"))
            wetp <- lapply(etp, `[`, c(info$dim_nms[1], "entropy_within"))
            etp <- lapply(etp, `[`, c(info$dim_nms[1], "entropy"))
        }
    } else {
        etp <- apply(occt, 1, entropy, log_base)
        if(!is.null(Theil_decomp)) {
            g <- apply(occt, 1, by, Theil_decomp, sum)
            g <- aperm(g, c(2, 1))
            betp <- apply(g, 1, function(x) entropy(x, log_base))
            wetp <- etp - betp
        }

    }

    #-- Final steps
    etp <- wide_to_long(etp, info$n_dim, info$dim_nms[info$nd[[1]]], "entropy")
    if(!is.null(Theil_decomp)) {
        wetp <- wide_to_long(wetp, info$n_dim,
                             info$dim_nms[info$nd[[1]]], "entropy_within")
        betp <- wide_to_long(betp, info$n_dim,
                             info$dim_nms[info$nd[[1]]], "entropy_between")
        etp <- Reduce(merge, list(etp, wetp, betp))
    }
    if (info$n_dims == 3) {
        # TODO
        # This is a redundancy needed because of the merge
        etp <- etp[order(etp[, 2], etp[, 1]), ]
    }
    return(etp)
}
