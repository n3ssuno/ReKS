#' @name
#' shannonEvenness
#'
#' @title
#' Shannon Evenness
#'
#' @description
#' The function returns the Shannon information evenness of a given table.
#' Optionally, it can return also the values for the evenness decomposed in
#' the two constituents of within- and between-groups evenness.
#'
#' @details
#' The function returns the Shannon information evenness of a given table.
#' It computes internally the frequency of each observation and, if not
#' specified differently, uses "bits" as units of measure.
#' The function can compute either the evenness of a cross-section, or of
#' a panel.
#' If a vector of groups is provided, the function returns also the values
#' for the evenness decomposed in the two constituents of within- and
#' between-groups evenness.
#'
#' @references
#' Pielou (1969) \emph{An Introduction to Mathematical Ecology}, Wiley.
#'
#' Stirling (2007) ``A General Framework for Analysing Diversity in Science,
#' Technology and Society'', \emph{Journal of The Royal Society Interface}, 4,
#' 707--719.
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
#' @return data.table with the (total) evenness and (eventually) its
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
#' EVS <- shannonEvenness(octab,
#'                       log_base = exp(1),
#'                       Theil_decomp = grps)
#'
#' @export

shannonEvenness <- function(occt,
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
    evenness <- function(d, log_base) {
        d <- d[which(d > 0)]
        f <- d / sum(d)
        e <- -sum(f * log(f, log_base))
        e[!is.finite(e)] <- 0
        e <- e / log(length(f), log_base)
        return(e)
    }

    #-- Main Function
    if (info$n_dims == 3) {
        evs <- apply(occt, 3, shannonEvenness, log_base, Theil_decomp)
        if(!is.null(Theil_decomp)) {
            bevs <- lapply(evs, `[`, c(info$dim_nms[1], "evenness_between"))
            wevs <- lapply(evs, `[`, c(info$dim_nms[1], "evenness_within"))
            evs <- lapply(evs, `[`, c(info$dim_nms[1], "evenness"))
        }
    } else {
        evs <- apply(occt, 1, evenness, log_base)
        if(!is.null(Theil_decomp)) {
            etp <- apply(occt, 1, function(x) entropy(x, log_base))
            g <- apply(occt, 1, by, Theil_decomp, sum)
            g <- aperm(g, c(2, 1))
            betp <- apply(g, 1, function(x) entropy(x, log_base))
            wetp <- etp - betp
            wevsmax <- apply(occt, info$nd[[1]],
                             function(x) log(length(x[which(x > 0)]), log_base))
            wevs <- wetp / wevsmax
            bevs <- apply(g, info$nd[[1]],
                          function(x) evenness(x, log_base))
        }
        # evs <- apply(occt, nd, function(x) entropy(x, log_base,
        #                                            evenness = TRUE))
        # if(!is.null(Theil_decomp)) {
        #     g <- apply(occt, nd, by, Theil_decomp, sum)
        #     g <- aperm(g, ndr)
        #     bevs <- apply(g, nd, function(x) entropy(x, log_base,
        #                                              evenness = TRUE))
        #     # wevs <- apply(octab, nd, by, grps, entropy,
        #     #               log_base, evenness = TRUE)
        #     # if(n_dims == 3) {
        #     #     wevs <- apply(wevs, 3, colSums)
        #     # } else {
        #     #     wevs <- colSums(wevs)
        #     # }
        #     wevsmax <- apply(g, nd, function(x) length(x[which(x > 0)]))
        #     #merge <- apply(wetp, 1, function(x) mapply(`/`, x, wevsmax))
        #     wevs <- Reduce(`/`, list(wetp, wevsmax))
        # }
    }

    #-- Final steps
    evs <- wide_to_long(evs, info$n_dim, info$dim_nms[info$nd[[1]]], "evenness")
    if(!is.null(Theil_decomp)) {
        wevs <- wide_to_long(wevs, info$n_dim,
                             info$dim_nms[info$nd[[1]]], "evenness_within")
        bevs <- wide_to_long(bevs, info$n_dim,
                             info$dim_nms[info$nd[[1]]], "evenness_between")
        evs <- Reduce(merge, list(evs, wevs, bevs))
    }
    if (info$n_dims == 3) {
        # TODO
        # This is a redundancy needed because of the merge
        evs <- evs[order(evs[, 2], evs[, 1]), ]
    }
    return(evs)
}
