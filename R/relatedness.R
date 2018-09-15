
# library(data.table)
# load("../OECD_REGPAT_2018/REGPAT.EU28.NUTS2.RData")
# d <- unique(regpat.eu28.nuts2[Prio_Year %in% 1996:2000,
#                               .(App_nbr, Prio_Year, IPC)])
# d <- d[, .(App_nbr, IPC = substr(IPC, 1, 7))]
# d <- xtabs(~ App_nbr + IPC, d, sparse = T)

# library(data.table)
# load("../OECD_REGPAT_2018/REGPAT.EU28.NUTS2.RData")
# d <- unique(regpat.eu28.nuts2[, .(App_nbr, IPC)])
# d <- d[, .(App_nbr, IPC = substr(IPC, 1, 7))]
# d <- xtabs(~ App_nbr + IPC, d, sparse = T)

relatedness <- function(adj_mtx, output_statistic = "t",
                        is_binary = NULL,
                        fixedmar = "both", seed = Sys.time(), nSim = 1000) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop(paste0("Package \"Matrix\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    # Preliminary operations

    adj_mtx <- as(adj_mtx, "Matrix")

    geo_dim <- attr(adj_mtx, 'geo_dim')
    kng_dim <- attr(adj_mtx, 'kng_dim')

    # t-stat

    relatedness_t <- function(...) {

        rnms <- rownames(adj_mtx)
        cnms <- colnames(adj_mtx)

        # adj_mtx <- as(adj_mtx, "ngCMatrix")
        adj_mtx[which(adj_mtx > 0)] <- 1
        if (any(adj_mtx@x != 1))
            stop(paste("It is not possible to transform the matrix into ",
                       "a binary (0/1) one"))

        Nr <- nrow(adj_mtx)

        J <- Matrix::crossprod(adj_mtx)
        Matrix::diag(J) <- 0

        n <- Matrix::colSums(adj_mtx)
        mu <- Matrix::tcrossprod(n)
        mu <- mu / Nr

        s2 <- Matrix::tcrossprod((1 - (n / Nr)), ((Nr - n) / (Nr - 1)))
        s2 <- mu * s2

        t <- (J - mu) / sqrt(s2)

        Matrix::diag(t) <- 0

        rownames(t) <- colnames(t) <- cnms

        return(t)
    }

    # p-value

    relatedness_p <- function(...) {
        if (!requireNamespace("vegan", quietly = TRUE)) {
            stop(paste0("Package \"vegan\" needed for this function to work. ",
                        "Please install it."), call. = FALSE)
        }

        rnms <- rownames(adj_mtx)
        cnms <- colnames(adj_mtx)

        isBinary <- ifelse((is.null(is_binary) &&
                                all(adj_mtx@x %in% c(0, 1, FALSE, TRUE))) ||
                               isTRUE(is_binary), TRUE, FALSE)

        if (isBinary)
            adj_mtx <- as(adj_mtx, "ngCMatrix")

        set.seed(seed)
        adj_mtx_null_models <- vegan::permatswap(adj_mtx,
                                                 fixedmar = fixedmar,
                                                 mtype = ifelse(isBinary,
                                                                "prab",
                                                                "count"),
                                                 times = nSim)

        J_hat <- Matrix::crossprod(adj_mtx)
        J <- lapply(adj_mtx_null_models$perm,
                    function(m) Matrix::crossprod(as(m, class(adj_mtx)[[1]])))
        p <- lapply(J, function(m) J_hat >= m)
        p <- Reduce("+", p)
        p <- p / nSim

        Matrix::diag(p) <- 0

        pPlus <- pmax(as.vector(2 * p - 1), rep(0, nrow(p) * ncol(p)))
        pPlus <- Matrix(pPlus, nrow = nrow(p), ncol = ncol(p))
        pMinus <- pmin(as.vector(2 * p - 1), rep(0, nrow(p) * ncol(p))) * (-1)
        pMinus <- Matrix(pMinus, nrow = nrow(p), ncol = ncol(p))

        rownames(p) <- colnames(p) <- cnms
        rownames(pPlus) <- colnames(pPlus) <- cnms
        rownames(pMinus) <- colnames(pMinus) <- cnms

        return(list(p = p, pPlus = pPlus, pMinus = pMinus))
    }

    R <- switch(output_statistic,
                t = relatedness_t(),
                p = relatedness_p(),
                stop('\"output_statistic\" can be one of \"t\", or \"p\"'))

    # class(R) <- c("reks_relatedness", class(R))
    attr(R, output_statistic) <- output_statistic
    attr(R, 'geo_dim') <- geo_dim
    attr(R, 'kng_dim') <- kng_dim

    return(R)
}
