.get_biadj_matrix <- function(data, geo_dim, kng_dim, kng_nbr,
                               binary_mode = "none") {
    # It returns a binary biadjacency matrix
    # It is possible to use the RTAs as a criterion of binarisation or simply
    #  to assume that each observation produces a 1 in the matrix, regardless of
    #  its "value"

    # Possible binary_mode options:
    # - none: weighted matrix (not binary);
    # - simple: if there is at least a patent in the geo. area in the
    #    particular tech. class it is 1, and 0 otherwise;
    # - RTA or RCA: Balassa method;
    # - higher_quartiles: similar to 'simple', but it exludes the first quartile
    #    (it seems the best choice);
    # - higher_quartiles_kng: similar to 'higher_quartiles', but the quartiles
    #    are computed for each knowledge class, assuming that some are more
    #    "ubiquitous" than others;
    # - higher_deciles_kng: similar to 'higher_quartiles_kng', but it excludes
    #    the lower 10% rather than the lower 25%.

    # Preliminary controls ---------------

    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop(paste0("Package \"reshape2\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    data <- as.data.frame(data)

    # Main function ---------------
    # if (is.null(kng_nbr)) {
    #     BM <- unique(data[, c(geo_dim, kng_dim)])
    #     BM$C <- 1
    #     BM <- reshape2::dcast(BM,
    #                           formula(paste(geo_dim, "~", kng_dim)),
    #                           value.var = "C")
    # } else {
    BM <- unique(data[, c(geo_dim, kng_nbr, kng_dim)])
    BM <- reshape2::dcast(BM,
                          formula(paste(geo_dim, "~", kng_dim)),
                          value.var = kng_nbr,
                          fun.aggregate = sum)
    # }

    BM[is.na(BM)] <- 0
    gnames <- BM[, 1]
    BM <- as.matrix(BM[, -1])
    rownames(BM) <- gnames

    if (binary_mode == 'RTA' | binary_mode == 'RCA') {
        BM <- .get_rta(BM, binary = TRUE)
    }
    if (binary_mode == 'simple') {
        BM[BM > 0] <- 1
    }
    if (binary_mode == "higher_quartiles") {
        fq <- quantile(BM[BM > 0])[2]
        BM[BM <  fq] <- 0
        BM[BM >= fq] <- 1
    }
    if (binary_mode == "higher_quartiles_kng") {
        fqs <- apply(BM, 2, function(col) quantile(col[col > 0])[2])
        knames <- colnames(BM)
        BM <- sapply(1:ncol(BM), function(col) {
            bm <- BM[, col]
            bm[bm <  fqs[col]] <- 0
            bm[bm >= fqs[col]] <- 1
            return(bm)
        })
        colnames(BM) <- knames
    }
    if (binary_mode == "higher_deciles_kng") {
        fqs <- apply(BM, 2,
                     function(col) quantile(col[col > 0], seq(0, 1, .1))[2])
        knames <- colnames(BM)
        BM <- sapply(1:ncol(BM), function(col) {
            bm <- BM[, col]
            bm[bm <  fqs[col]] <- 0
            bm[bm >= fqs[col]] <- 1
            return(bm)
        })
        colnames(BM) <- knames
    }

    class(BM) <- c('reks_biadj_matrix', 'matrix')
    attr(BM, "geo_dim") <- geo_dim
    attr(BM, "kng_dim") <- kng_dim
    attr(BM, "binary")  <- ifelse(binary_mode == "none", FALSE, TRUE)
    if (binary_mode == 'RTA') {
        attr(BM, "binary_mode") <- 'RTA'
    }
    if (binary_mode == 'RCA') {
        attr(BM, "binary_mode") <- 'RCA'
    }
    if (binary_mode == 'simple') {
        attr(BM, "binary_mode") <- 'simple'
    }
    if (binary_mode == "higher_quartiles") {
        attr(BM, "binary_mode") <- 'higher_quartiles'
    }
    if (binary_mode == "higher_quartiles_kng") {
        attr(BM, "binary_mode") <- 'higher_quartiles_kng'
    }
    if (binary_mode == "higher_deciles_kng") {
        attr(BM, "binary_mode") <- 'higher_quartiles_kng'
    }

    return(BM)
}
