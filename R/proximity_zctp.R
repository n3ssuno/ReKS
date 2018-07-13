#' @name
#' proximity_zctp
#'
#' @title
#' Proximity a là Zaccaria, Cristelli, Tacchella, Pietronero
#'
#' @description
#' The function computes the so called Proximity, as proposed by Zaccaria,
#' Cristelli, Tacchella, and Pietronero in 2014 on \emph{PLoS ONE}.
#'
#' @details
#'
#'
#' @references
#' Zaccaria, Cristelli, Tacchella and Pietronero (2014) ``How the Taxonomy of
#' Products Drives the Economic Development of Countries'', \emph{PLoS ONE},
#' 0113770, 1--17.
#'
#' Zhou, Ren, Medo and Zhang (2007) ``Bipartite Network Projection and Personal
#' Recommendation'', \emph{Physical Review E}, 76, 046115.
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
#' @return A data.frame with the Coherence Index of each geographical area.
#'
#' @examples
#' RCI <- coherence(data = df, geo_dim = NUTS2,
#' kng_dim = IPC.3dig, kng_nbr = N.patents)

proximity_zctp <- function (data, geo_dim, kng_dim, kng_nbr,
                            binary_mode = "RTA", names_as_strings = FALSE) {

    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop(paste0("Package \"reshape2\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    data <- as.data.frame(data)

    if (names_as_strings == FALSE) {
        geo_dim <- deparse(substitute(geo_dim))
        kng_dim <- deparse(substitute(kng_dim))
        kng_nbr <- deparse(substitute(kng_nbr))
    }

    data <- data[complete.cases(data[, c(geo_dim, kng_dim, kng_nbr)]), ]

    mm <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr,
                            binary_mode)

    if (any(rowSums(mm) == 0)) {
        mm <- mm[-which(rowSums(mm) == 0), ]
    }
    if (any(colSums(mm) == 0)) {
        mm <- mm[, -which(colSums(mm) == 0)]
    }

    du <- .get_du(mm)

    mm_proj <- t(mm) %*% (mm / du$diversification)

    max_u <- apply(expand.grid(du$ubiquity, du$ubiquity), 1, max)
    max_u <- matrix(max_u, nrow = length(du$ubiquity))

    m_out <- mm_proj * (1/max_u)

    if (!isSymmetric.matrix(m_out)) {
        stop(paste("Something went wrong: \n",
                   "the matrix produced is not simmetric as it should be!"))
    }
    if (!all(m_out>=0 & m_out<=1, na.rm = T)) {
        stop(paste("Something went wrong: \n",
                   "at least one value is not a probability ",
                   "(i.e., x !in [0,1])!"))
    }

    diag(m_out) <- NA
    #m_out[upper.tri(m_out)] <- NA

    B <- reshape2::melt(m_out, na.rm = TRUE,
                        varnames = c(paste0(kng_dim, "_1"),
                                     paste0(kng_dim, "_2")),
                        value.name = "B")

    class(B) <- c('data.frame', 'reks_proximity_zctp')
    attr(B, 'diversification') <- du$diversification
    attr(B, 'ubiquity') <- du$ubiquity
    if (binary_mode == 'RTA') {
        attr(B, "binary_mode") <- 'RTA'
    }
    if (binary_mode == 'RCA') {
        attr(B, "binary_mode") <- 'RCA'
    }
    if (binary_mode == 'simple') {
        attr(B, "binary_mode") <- 'simple'
    }
    if (binary_mode == "higher_quartiles") {
        attr(B, "binary_mode") <- 'higher_quartiles'
    }
    if (binary_mode == "higher_quartiles_kng") {
        attr(B, "binary_mode") <- 'higher_quartiles_kng'
    }
    attr(B, "geo_dim") <- geo_dim
    attr(B, "kng_dim_1") <- paste0(kng_dim, "_1")
    attr(B, "kng_dim_2") <- paste0(kng_dim, "_2")
    attr(B, "kng_nbr") <- kng_nbr
    attr(B, "prx") <- "B"
    attr(B, "directed") <- TRUE

    return(B)

}
