.get_du <- function(data) {
    du <- list()
    du$diversification <- rowSums(mm)
    du$ubiquity  <- colSums(mm)
    return(du)
}
