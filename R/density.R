
local_density <- function(occt) {
    #-- Preliminary steps and checks
    info <- data_info(occt)
    occt <- remove_zeros(occt)

    #-- Main function
    if (info$n_dims == 3) {
        w <- apply(occt, 3, local_density)
    } else {
        Phi <- proximity(occt)
        den <- Matrix::rowSums(Phi, na.rm = TRUE)
        Mcp <- rta(occt, binary = TRUE)
        w <- Matrix::tcrossprod(Mcp, Phi / den)
    }

    return(w)
}

