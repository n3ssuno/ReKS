#TODO
#doesn't work for reks_entropy class because there's more than one measure

.get_cor <- function(data, method = 'spearman') {
    # TODO which measures
    if (!any(class(data) %in% c("reks_coherence",
                                "reks_hh_complexity",
                                "reks_fitness_tccgp"))) {
        # TODO error message
        stop('')
    }

    data <- as.data.frame(data)

    #TODO add in the measures
    geo_dim <- attr(data, 'geo_dim')
    time_dim <- attr(data, 'time_dim')
    measure <- attr(data, 'measure')

    time_span <- unique(data[, time_dim])

    data <- dcast(data,
                  paste(geo_dim, "~", time_dim),
                  value.var = measure)
    rownames(data) <- data[, 1]
    data <- data[, -1]

    COR <- lapply(1:(length(time_span) - 1), function(t) {
        c(firstYear = time_span[t],
          sectonYear = time_span[t + 1],
          cor = cor(data[, t], data[, t + 1],
                    use = 'pairwise.complete.obs', method = method))
    })
    COR <- do.call('rbind', COR)
    COR <- as.data.frame(COR)

    class(COR) <- c('reks_correlation_coefficient', 'data.frame')
    attr(cor_coef, 'measure') <- measure
    attr(cor_coef, 'method') <- method

    return(COR)
}

plot.reks_correlation_coefficient <- function(cor_coef) {
    measure <- attr(cor_coef, 'measure')
    plot(cor_coef[, "firstYear"], cor_coef[, "cor"],
         xlab = "Year", ylab = "Correlation",
         main = paste(measure,
                      "- Spearman correlation between consecutive years"))
}

plot_correlation_coefficient <- function(data) {
    COR <- .get_cor(data)
    plot(COR)
}
