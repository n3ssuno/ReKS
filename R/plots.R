
plot_biadjacency_matrix <- function (occt, order = "DU", ...) {
    nms <- names(dimnames(occt))
    occt <- Matrix::Matrix(occt, sparse = TRUE)
    if (order == "DU") {
        diversity <- Matrix::rowSums(occt)
        ubiquity <- Matrix::colSums(occt)
        row_order <- order(diversity, decreasing = T)
        col_order <- order(ubiquity)
    }
    occt <- occt[row_order, col_order]
    Matrix::image(Matrix::t(occt),
                  xlab = nms[1],
                  ylab = nms[2],
                  ...)
}
