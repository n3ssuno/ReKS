# https://www.r-bloggers.com/creating-an-image-of-a-matrix-in-r-using-image/

plot_biadj_matrix <- function(data, geo_dim, kng_dim, kng_nbr,
                              binary_mode, order = "DU") {
    # Order can be
    # - "DU" = diversification - ubiquity
    # - "FC" = fitness - complexity
    # - "HH" = Hidalgo-Hausmann complexity
    # TODO
    # The last two haven't been implemented yet
    data <- as.data.frame(data)
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim <- deparse(substitute(kng_dim))
    kng_nbr <- deparse(substitute(kng_nbr))
    BM <- .get_biadj_matrix(data, geo_dim, kng_dim, kng_nbr, binary_mode)
    plot(BM, order)
}

plot.rks_biadj_matrix <- function(BM, order = "DU") {
    if (order == "DU") {
        du <- .get_du(BM)
        row_order <- order(du$diversification, decreasing = T)
        col_order <- order(du$ubiquity, decreasing = T)
    }
    bm <- BM[row_order, col_order]
    rotate <- function(x) t(apply(x, 2, rev))
    image(rotate(bm),
          col = c('white', 'black'),
          xlab = attr(BM, "kng_dim"), ylab = attr(BM, "geo_dim"),
          axes = FALSE)
    box()
}
