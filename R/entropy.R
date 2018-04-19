#' @name
#' entropy
#'
#' @title
#' Information Entropy
#'
#' @description
#' The function will return the frequency distribution of a given vector.
#' It is useful for entropy functions.
#'
#' @details
#'
#' @param data Frequency distribution of a vector of elements.
#' @return The information entropy level of the vector.
#' @examples
#' frequences <- get_freqs(c(1,3,5,2,4,9,1,6,10))
#' etpy <- entropy(frequences)

entropy <- function(data) {
    freqs <- get_freqs(data)
    etp <- -sum(freqs * log2(freqs))
    etp[!is.finite(etp)] <- 0
    return(etp)
}
