#' @name
#' relatedness
#'
#' @title
#' Relatedness
#'
#' @description
#' The function computes the so called Relatedness index.
#'
#' @details
#' The function computes the so called Relatedness index.
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
#' 1289-1302;
#'
#' @encoding UTF-8
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
#' @param output_statistic It can be "t" for t-statistic or "p" for p-value
#' (look at Bottazzi and Pirino 2010 for an explanation).
#' @param is_binary It can be NULL (default), TRUE or FALSE.
#' @param fixedmar It is useful to choose which constraint you want to impose
#' in the null model.
#' @param seed It is useful to choose the seed.
#' @param n_sim You can choose how many simulations you want to run to compute
#' the null model.
#' @param sparse It can be TRUE or FALSE
#'
#' @return A matrix of similarity between technological domains
#'
#' @export

# @importFrom methods as

# - merge the two functions
# - add other similarity measures (like Zhou et al 2007 and
#   Zaccaria et al. 2014)
# - improve the documentation

relatedness <- function(occt,
                        output_statistic = "t",
                        is_binary = NULL,
                        fixedmar = "both",
                        seed = Sys.time(),
                        n_sim = 1000,
                        sparse = FALSE) {

    #-- Preliminary steps and checks
    info <- data_info(occt)

    if (info$n_dims == 3) {
        occt <- apply(occt, c(1, 2), sum)
    }

    occt <- Matrix::Matrix(occt, sparse = sparse)

    # Internal functions ----
    # t-stat --------
    relatedness_t <- function(...) {
        if (all(occt > 0)) {
            stop(paste("The adjacency matrix is complete.",
                       "It is not possible to compute the relatedness",
                       "between its nodes."))
        }
        occt[Matrix::which(occt > 0)] <- 1
        if (!all(occt@x %in% c(0, 1)))
            stop(paste("It is not possible to transform the matrix into ",
                       "a binary (0/1) one"))
        Nr <- nrow(occt)
        J <- Matrix::crossprod(occt)
        Matrix::diag(J) <- 0
        n <- Matrix::colSums(occt)
        mu <- Matrix::tcrossprod(n)
        mu <- mu / Nr
        s2 <- Matrix::tcrossprod(1 - (n / Nr), (Nr - n) / (Nr - 1))
        s2 <- mu * s2
        t <- (J - mu) / sqrt(s2)
        Matrix::diag(t) <- 0
        return(t)
    }
    # p-value --------
    relatedness_p <- function(...) {
        if (!requireNamespace("vegan", quietly = TRUE)) {
            stop(paste0("Package \"vegan\" needed for this function to work. ",
                        "Please install it."), call. = FALSE)
        }
        is_binary <- ifelse((is.null(is_binary) &&
                                all(occt@x %in% c(0, 1, FALSE, TRUE))) ||
                               isTRUE(is_binary), TRUE, FALSE)
        #if (is_binary)
        #    occt <- as(occt, "ngCMatrix")

        set.seed(seed)
        occt_null_models <- vegan::permatswap(occt,
                                              fixedmar = fixedmar,
                                              mtype = ifelse(is_binary,
                                                             "prab",
                                                             "count"),
                                              times = n_sim)

        J_hat <- Matrix::crossprod(occt)
        J <- lapply(occt_null_models$perm, function(m) {
            m <- Matrix::Matrix(m, sparse = sparse)
            m <- Matrix::crossprod(m)
            return(m)
            })
        p <- lapply(J, function(m) J_hat >= m)
        p <- Reduce("+", p)
        p <- p / n_sim

        Matrix::diag(p) <- 0

        p_plus <- pmax(as.vector(2 * p - 1), rep(0, nrow(p) * ncol(p)))
        p_plus <- Matrix::Matrix(p_plus, nrow = nrow(p), ncol = ncol(p))
        p_minus <- pmin(as.vector(2 * p - 1), rep(0, nrow(p) * ncol(p))) * (-1)
        p_minus <- Matrix::Matrix(p_minus, nrow = nrow(p), ncol = ncol(p))

        names(attr(p, "Dimnames")) <- info$dim_nms[info$nd[[2]]]
        names(attr(p_plus, "Dimnames")) <- info$dim_nms[info$nd[[2]]]
        names(attr(p_minus, "Dimnames")) <- info$dim_nms[info$nd[[2]]]

        #rownames(p) <- colnames(p) <- cnms
        #rownames(p_plus) <- colnames(p_plus) <- cnms
        #rownames(p_minus) <- colnames(p_minus) <- cnms

        attr(p, "p_plus") <- p_plus
        attr(p, "p_minus") <- p_minus

        return(p)
    }

    R <- switch(output_statistic,
                t = relatedness_t(),
                p = relatedness_p(),
                stop('\"output_statistic\" can be one of \"t\", or \"p\"'))

    # class(R) <- c("reks_relatedness", class(R))
    #attr(R, output_statistic) <- output_statistic
    #attr(R, 'geo_dim') <- geo_dim
    #attr(R, 'kng_dim') <- kng_dim

    return(R)
}

#' @name
#' proximity
#'
#' @title
#' Proximity
#'
#' @description
#' The function computes the so called Proximity index.
#'
#' @details
#' The function computes the so called Proximity index.
#'
#' @references
#' Hidalgo, Klinger, Barab{\'a}si, and Hausmann (2007) ``The Product Space
#' Conditions the Development of Nations'', \emph{Science}.
#'
#' @encoding UTF-8
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
#'
#' @return A matrix of similarity between technological domains

# - Works only for 2D objects
# - improve documentation

proximity <- function(occt) {

    occt <- remove_zeros(occt)

    Mcp <- rta(occt, binary = TRUE)
    Phi <- Matrix::crossprod(Mcp)
    ubiquity <- Matrix::colSums(Mcp)
    max_u <- apply(expand.grid(ubiquity, ubiquity), 1, max)
    max_u <- matrix(max_u, nrow = length(ubiquity))
    Phi <- Phi / max_u

    #Phi <- (1 + cor(as.matrix(RTA))) / 2

    return(Phi)
}

#' @name
#' association_strength
#'
#' @title
#' Association Strength

# - Works only for 2D objects
# - improve documentation

association_strength <- function(occt, SCALE = NULL) {

    occt <- remove_zeros(occt)

    Phi <- Matrix::crossprod(occt)
    Phi <- rta(Phi)

    if (!is.null(SCALE)) {
        Phi <- SCALE(Phi)
    }

    return(Phi)
}
