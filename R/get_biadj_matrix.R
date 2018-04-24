.get_biadj_matrix <- function(data, geo_dim, kng_dim,
                               RTA = FALSE, kng_nbr = NULL) {
    # It returns a binary biadjacency matrix
    # It is possible to use the RTAs as a criterion of binarisation or simply
    #  to assume that each observation produces a 1 in the matrix, regardless of
    #  its "value"

    # Preliminary controls ---------------
    if (isTRUE(RTA) & (is.null(kng_nbr) || kng_nbr == 'NULL')) {
        stop(paste('You must provide a kng_nbr column name in the options ',
                   'if you want to compute the Relative Technological ',
                   'Advantages of each geographical area.'))
    }
    # if (!is.null(kng_nbr) & !isTRUE(RTA)) {
    #     stop(paste('By now, it is possible to use only the Relative ',
    #                'Technological Advantages as a way to have binary matrices ',
    #                'from non-binary ones. Please, choose RTA = TRUE in the ',
    #                'options.'))
    # }
    data <- as.data.frame(data)

    # Main function ---------------
    if (is.null(kng_nbr)) {
        BM <- unique(data[, c(geo_dim, kng_dim)])
        BM$C <- 1
        BM <- reshape2::dcast(BM,
                              formula(paste(geo_dim, "~", kng_dim)),
                              value.var = "C")
    } else {
        BM <- unique(data[, c(geo_dim, kng_nbr, kng_dim)])
        BM <- reshape2::dcast(BM,
                              formula(paste(geo_dim, "~", kng_dim)),
                              value.var = kng_nbr,
                              fun.aggregate = sum)
    }

    BM[is.na(BM)] <- 0
    gnames <- BM[, 1]
    BM <- as.matrix(BM[, -1])
    rownames(BM) <- gnames

    if (isTRUE(RTA)) {
        BM <- .get_rta(BM, binary = TRUE)
    }

    class(BM) <- c('rks_biadj_matrix', 'matrix')
    attr(BM, "geo_dim") <- geo_dim
    attr(BM, "kng_dim") <- kng_dim
    attr(BM, "binary")  <- TRUE
    attr(BM, "RTA")     <- TRUE

    return(BM)
}
