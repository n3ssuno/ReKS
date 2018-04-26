.get_du <- function(biadj_matrix) {
    du <- list()
    du$diversification <- rowSums(biadj_matrix)
    du$ubiquity  <- colSums(biadj_matrix)
    return(du)
}
