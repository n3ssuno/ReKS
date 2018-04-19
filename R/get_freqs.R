#' @name
#' get_freqs
#'
#' @title
#' Get the frequency distribution of a given vector
#'
#' @description
#' The function will return the frequency distribution of a given vector.
#' It is useful for entropy functions.
#'
#' @details
#'
#' @param data Vector of numbers. If it sums to 1, it assumes they are already
#' frequencies.
#' @return A vector of numbers that sums to 1.
#' @examples
#' frequences <- get_freqs(values)
#' frequences <- get_freqs(table(values))

get_freqs <- function(data) {
    if (sum(as.numeric(data))==1) {
        warning('I assume you provided me a list of relative frequences')
        return(data)
    } else {
        warning(paste('I assume you provided me a list of absolute values.\n',
                      'I internally transformed them in relative frequencies.\n',
                      'Otherwise check in the original data',
                      ' why their sum is not 1.'))
        freqs <- data/sum(as.numeric(data))
    }
    return(freqs)
}
