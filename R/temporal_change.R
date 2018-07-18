#' @title temporal_change
#' @param data Any onject produced by a ReKS "relatedness measure" function
#' @param delta A number or a vector of years (it will return the growth rate of
#' the measure chosen between a year and the same measure delta years later)
#' @return A data.frame (or a list of data.frames - with a warning in this case)
#' in which for each delta it has been computed the growth rate between the
#' measure in a given year and the same measure delta years later.

temporal_change <- function(data, delta = 1) {
    data <- as.data.frame(data)

    time_dim <- attr(data, 'time_dim')
    geo_dim <- attr(data, 'geo_dim')
    measure <- attr(data, 'measure')

    if (class(data[, time_dim]) == "factor") {
        tl <- levels(data[, time_dim])
        tn <- data[, time_dim]
        data[, time_dim] <- as.numeric(tl)[tn]
    }

    time_span <- unique(data[, time_dim])

    if (max(delta) >= length(time_span)) {
        stop(paste0("You cannot use a delta greater than ",
                    "the number of years considered"))
    }

    DDsM <- lapply(measure, function(m) {
        DDs <- lapply(delta, function(d) {
            DD <- lapply(time_span[1:(length(time_span) - d)], function(t) {
                t1 <- t
                t2 <- t + d
                firstYear <- data[which(data[, time_dim] == t1),
                                  c(geo_dim, m)]
                names(firstYear) <- c(geo_dim, t1)
                secndYear <- data[which(data[, time_dim] == t2),
                                  c(geo_dim, m)]
                names(secndYear) <- c(geo_dim, t2)
                dta <- merge.data.frame(firstYear, secndYear)
                t1n <- as.character(t1)
                t2n <- as.character(t2)
                D <- dta[, t2n] - dta[, t1n]
                D <- D / dta[, t1n]
                return(cbind.data.frame(dta[, geo_dim], t1n,
                                        #t2n,
                                        D))
            })
            DD <- do.call('rbind', DD)
            colnames(DD) <- c(geo_dim,
                              "firstYear",
                              #"secondYear",
                              paste0("Delta_", d))
            return(DD)
        })
        DDs <- plyr::join_all(DDs)

        return(DDs)
    })
    names(DDsM) <- measure

    class(DDsM) <- c("reks_temporal_change", class(DDsM))

    if (length(measure) == 1) {
        DDsM <- DDsM[[1]]
        attr(DDsM, 'measure') <- measure
        return(DDsM)
    } else {
        warning(paste0("It has been returned a list of objects: ",
                       "one for each measure considered"))
        return(DDsM)
    }
}

print.reks_temporal_change <- function(TC) {
    if (any(class(TC) == "list")) {
        par(mfrow = c(ceiling((length(TC) + 1)/2), 2))
        nd <- ncol(TC[[1]])
        measure <- names(TC)
        cols <- rainbow(nd - 2)
        for (m in measure) {
            plot(density(TC[[m]][, "Delta_1"], na.rm = TRUE),
                 main = m, col = cols[1])
            if (nd >= 4) {
                for (d in (4:nd) - 2) {
                    lines(density(TC[[m]][, paste0("Delta_", d)], na.rm = TRUE),
                          main = m, col = cols[d])
                }
            }
        }
        plot.new()
        legend("center",
               legend = paste0("Delta_", 1:(nd - 2)),
               col = cols, lty = 1)
    } else {
        par(mfrow = c(1, 1))
        nd <- ncol(TC)
        measure <- attr(TC, "measure")
        cols <- rainbow(nd - 2)
        for (m in measure) {
            plot(density(TC[, "Delta_1"], na.rm = TRUE),
                 main = m, col = cols[1])
            if (nd >= 4) {
                for (d in (4:nd) - 2) {
                    lines(density(TC[, paste0("Delta_", d)], na.rm = TRUE),
                          main = m, col = cols[d])
                }
            }
        }
        legend("topright",
               legend = paste0("Delta_", 1:(nd - 2)),
               col = cols, lty = 1)
    }

}
