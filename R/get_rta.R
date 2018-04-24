.get_rta <- function(data, binary = FALSE) {
    # It returns the Relative Technological Advantages (i.e., Balassa RCAs) of
    #  each row of a matrix. The name comes from Antonelli et al. (2017).
    # If binary == TRUE, it returns, instead of the RTAs, a binary matrix
    #  in which each cell has value 1 if its RTA is >= 1, and 0 otherwise.
    # See:
    # - Balassa (1965) "Trade Liberalisation and Revealed Comparative
    #    Advantage", The Manchester School, 33, 99–123;
    # - Antonelli, Crespi, Mongeau Ospina and Scellato (2017) "Knowledge
    #    Composition, Jacobs Externalities and Innovation Performance in
    #    European Regions", Regional Studies, 51, 1708-1720.

    data <- try(as.data.frame(data))

    RTA <- t(t(data / rowSums(data)) / (colSums(data) / sum(data)))

    if (isTRUE(binary)) {
        RTA <- ifelse(RTA >= 1, 1, 0)
    }

    return(RTA)
}
