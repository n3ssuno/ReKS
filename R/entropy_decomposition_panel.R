#' @name
#' entropy_decomposition_panel
#'
#' @title
#' Entropy Decomposition Theorem for panel data
#'
#' @description
#' The function will return the so called Related/Unerlated Variety vales for
#' each region in each year.
#'
#' @details
#' The function will return the so called entropy decomposition theorem
#' so that, besides the overall entropy level, also the value of the between
#' and within group entropy is provided. Internally the function computes and
#' temorarily stores also other values. However they are not provided in the
#' oputput DB.
#'
#' @param data It is expected to be a data.frame, data.table or an array of
#' numbers. This is the only mandatory input, all others being optional.
#' @param kng_nbr It is expected to be the name of the (numeric) column of the DB
#' in which there is either the absolute number of patents (firms) of a given
#' class (industry), or the corresponding relative frequences.
#' @param geo_dim It is expected to be a (character) column of the DB that
#' identifies the first level of aggregation of the data. E.g., can be a
#' geographical area (county, region, nation) or can refer to the temporal
#' dimension (year). The name of the column in the output DB can be chosen
#' through the geo_dim_label option, otherwise a generic label will be use.
#' @param time_dim It is expected to be a (character) column of the DB that
#' identifies the second level of aggregation of the data. It is take for
#' granted that you cannot use it without geo_dim.
#' @param kng_dim_upper It is expected to be an integer that denotes the number
#' of digits considered as the level of aggregation to compute the entropy
#' decomposition (that is, these are the groups to which the word 'within
#' groups' and 'between groups' refer)
#' @return the value of the entropy index if you provide only an array of
#' values. Otherwise, a data.table with the overall entropy index, and its
#' decomposition in between- and within-group components.
#' @examples
#' Variety(data = data, kng_nbr = N, geo_dim = region,
#' time_dim = year, kng_dim_upper = IPC.1,
#' geo_dim_label = "region", time_dim_label = "year")
#' Variety(data = data$N)
#' Variety(data = data, kng_nbr = freq, kng_dim_upper = IPC.3)

entropy_decomposition_panel <- function(data, kng_nbr, kng_dim_upper,
                                        geo_dim, time_dim) {

    # Preliminary transformations and checks -------------

    data <- as.data.frame(data)
    if (!all(complete.cases(data))) {
        warning(paste('There is some non complete row in the database.\n',
                      'I cannot guarrenty you about the results.'))
    }

    kng_nbr <- deparse(substitute(kng_nbr))
    geo_dim <- deparse(substitute(geo_dim))
    kng_dim_upper <- deparse(substitute(kng_dim_upper))
    time_dim <- deparse(substitute(time_dim))

    # Decomposed entropy -------------

    dd <- split(data[, kng_nbr],
                list(data[, time_dim], data[, geo_dim]))
    ddnt <- sapply(names(dd), function(s) strsplit(s, "[.]")[[1]][1])
    ddng <- sapply(names(dd), function(s) strsplit(s, "[.]")[[1]][2])

    obs_list <- 1:length(dd)
    entropy_total <- sapply(obs_list,
                              function(x) RKS::entropy(dd[[x]]))
    entropy_total <- cbind.data.frame(ddnt, ddng, entropy_total)
    colnames(entropy_total) <- c(time_dim, geo_dim, "entropy.total")

    grps <- split(data[, kng_dim_upper],
                  list(data[, time_dim], data[, geo_dim]))
    entropy_decomposed <- sapply(obs_list,
                                 function(x)
                                     RKS::entropy_decomposition(dd[[x]],
                                                                grps[[x]]))
    entropy_decomposed <- matrix(unlist(entropy_decomposed[1:2,]),
                                 ncol = 2, byrow = T)
    entropy_decomposed <- cbind.data.frame(ddnt, ddng, entropy_decomposed)
    colnames(entropy_decomposed) <- c(time_dim, geo_dim,
                                      "entropy.between", "entropy.within")

    entropy <- merge(entropy_total, entropy_decomposed)

    return(entropy)
}
