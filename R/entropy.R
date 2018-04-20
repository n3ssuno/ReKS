#' @name
#' entropy
#'
#' @title
#' Information Entropy
#'
#' @description
#' The function will return the information entropy of a given vector that
#' reports the absolute frequency of each of the possible types/groups of
#' observations of the database.
#'
#' @details The function will return the information entropy of a given vector
#' that reports the absolute frequency of each of the possible types/groups of
#' observations of the database.
#' It is also used internally by the entropy decomposition functions.
#' See:
#' \itemize{
#' \item{Shannon (1948) "A Mathematical Theory of Communication", \emph{Bell
#' System Technical Journal}, 27, 379--423;}
#' \item{Theil (1967) \emph{Economics and Information Theory}, North-Holland;}
#' \item{Theil (1972) \emph{Statistical Decomposition Analysis}, North-Holland;}
#' \item{Frenken (2007) "Entropy statistics and information theory", in
#' Hanusch and Pyka (Eds.) \emph{Elgar Companion to Neo-Schumpeterian
#' Economics}, Edward Elgar.}
#' }
#'
#' @param data Numerical vector whos elements are the numerosity of each
#' of the types in which the observations can be categorised (absolute
#' frequences).
#' @return The information entropy level of the vector.
#' @examples
#' etpy <- entropy(c(1,3,5,2,4,9,1,6,10))

entropy <- function(data) {
    freqs <- .get_freqs(data)
    etp <- -sum(freqs * log2(freqs))
    etp[!is.finite(etp)] <- 0
    return(etp)
}
