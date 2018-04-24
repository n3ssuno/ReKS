plot_biadj_matrix <- function(data, geo_dim, kng_dim, order = "DU") {
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
    BM <- .get_biadj_matrix(data, geo_dim, kng_dim)
    plot(BM)
}

plot.rks_biadj_matrix <- function(BM) {
    bm <- BM[order(rowSums(BM)), # use diversity
              order(colSums(BM), decreasing = T)] # use ubiquity
    image(bm,
          col = c('black', 'white'),
          xlab = attr(BM, "kng_dim"), ylab = attr(BM, "geo_dim"),
          axes = FALSE)
    box()
}
