# possible to use also fitness
# merge in a unique function coi_cog()
# recursive as the others and remove "long" option

#' @export

coi <- function(occt, scale = TRUE, long = TRUE) {
    # coi: Complexity Outlook Index

    info <- data_info(occt)

    if (info$n_dims == 3) {
        COI <- apply(occt, info$n_dims, coi, long = FALSE)
    } else {
        Mcp <- rta(occt, binary = TRUE)
        Mcp <- remove_zeros(Mcp)
        MM <- 1 - Mcp
        #prox_m <- proximity(occt)
        dens_m <- local_density(occt)
        tci <- complexity(as.table(occt), scale = FALSE, which = "tci")

        COI <- Matrix::rowSums(Matrix::t(Matrix::t(dens_m * MM) * tci$TCI))
        # check if it is correct comparing with py-complexity
        # if (scale == TRUE) {
        #     COI <- scale(COI)
        # }
    }

    if (long == TRUE) {
        COI <- wide_to_long(COI, info$n_dims, info$dim_nms[info$nd[[1]]], "coi")
    }

    return(COI)
}

#' @export

cog <- function(occt, scale = TRUE, long = TRUE) {
    # cog: Complexity Outlook Gain

    info <- data_info(occt)

    if (info$n_dims == 3) {
        COG <- apply(occt, info$n_dims, cog, long = FALSE)
        COG <- lapply(names(COG), function(x) cbind.data.frame(COG[[x]], x))
        COG <- Reduce(rbind, COG)
        COG <- COG[, c(1, 2, 4, 3)]
        colnames(COG) <- c(info$dim_nms, "cog")
    } else {
        Mcp <- rta(occt, binary = TRUE)
        Mcp <- remove_zeros(Mcp)
        MM <- 1 - Mcp
        prox_m <- proximity(occt)
        #dens_m <- local_density(occt)
        #rci <- complexity(as.table(occt), scale = FALSE, which = "rci")
        tci <- complexity(as.table(occt), scale = FALSE, which = "tci")

        COG <- MM *
            Matrix::tcrossprod(MM,
                               Matrix::t(prox_m *
                                             (tci$TCI /
                                                  Matrix::rowSums(prox_m))))
        # check if it is correct comparing with py-complexity
        # if (scale == TRUE) {
        #     COG <- COG / sd(rci)
        # }
        COG <- cbind.data.frame(expand.grid(dimnames(COG)),
                                cog = as.vector(COG))
    }

    return(COG)
}
