plot_maximum_spanning_tree <- function(phi_df,
                                       time_laps = NULL, min_thrashold = NULL,
                                       ...) {
    max_span_tree <- .maximum_spanning_tree(phi_df, time_laps, min_thrashold)
    plot(max_span_tree, ...)
}

#' @name plot.reks_maximum_spanning_tree
#'
#' @title Plot the Maximum Spanning Tree of a proximity matrix
#'
#' @description It computes the Maximum Spanning Tree for proximity matrices
#' of class \bold{reks_proximity_hkbh} or \bold{reks_proximity_zctp}, and plot
#' it.
#'
#' @details You can use the functions \code{proximity_hkbh()} or
#' \code{proximity_zctp()} (and the related ones for panel data) to get a
#' \code{phi_df} object of the proper class.
#'
#' @encoding UTF-8
#'
#' @param phi_df proximity matrix of type \bold{reks_proximity_hkbh} or
#' \bold{reks_proximity_zctp}
#' @param time_laps Year or vector of years to get a subset of the phi_df
#' (optional)
#' @param min_thrashold All the edges with a weight greater that this thrashold
#' will be added to the Maximum Spanning Tree
#' @return nothing
#'
#' @export

plot.reks_maximum_spanning_tree <- function(mst, ...) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(paste0("Package \"igraph\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    igraph::plot.igraph(mst, ...)
}
# vertex.label.cex = .5, vertex.size = .05


#' @name .maximum_spanning_tree
#'
#' @title Maximum Spanning Tree for proximity matrices
#'
#' @description It computes the Maximum Spanning Tree for proximity matrices
#' of class \bold{reks_proximity_hkbh} or \bold{reks_proximity_zctp}.
#'
#' @details You can use the functions \code{proximity_hkbh()} or
#' \code{proximity_zctp()} (and the related ones for panel data) to get a
#' \code{phi_df} object of the proper class.
#'
#' @encoding UTF-8
#'
#' @param phi_df proximity matrix of type \bold{reks_proximity_hkbh} or
#' \bold{reks_proximity_zctp}
#' @param time_laps Year or vector of years to get a subset of the phi_df
#' (optional)
#' @param min_thrashold All the edges with a weight greater that this thrashold
#' will be added to the Maximum Spanning Tree
#' @return Maximum Spanning Tree network (\bold{igraph} object)
#'
#' @keywords internal

.maximum_spanning_tree <- function(phi_df,
                                   time_laps = NULL, min_thrashold = NULL) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(paste0("Package \"igraph\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }

    # grph <- igraph::graph_from_adjacency_matrix(phi_mtx,
    #                                          mode = "lower",
    #                                          weighted = TRUE,
    #                                          diag = FALSE)

    proximity_classes <- c("reks_proximity_hkbh", "reks_proximity_zctp")
    if (!any(class(phi_df) %in% proximity_classes)) {
        stop(paste0("The object ", phi_df, " must belong to one of the ",
                    "following classes: ",
                    paste(proximity_classes, collapse = ", ")))
    }
    # geo_dim <- attr(phi_df, "geo_dim")
    kng_dim_1 <- attr(phi_df, "kng_dim_1")
    kng_dim_2 <- attr(phi_df, "kng_dim_2")
    # kng_nbr <- attr(phi_df, "kng_nbr")
    prx <- attr(phi_df, "prx")
    directed <- attr(phi_df, "directed")
    if (!is.null(attr(phi_df, "time_dim"))) {
        time_dim <- attr(phi_df, "time_dim")
        if (!all(time_laps %in% phi_df[, time_dim])) {
            stop(paste0("The time laps you have selected is not a subset of ",
                        "the time period covered by the dataset chosen"))
        }
        if (!is.null(time_laps)) {
            phi_df <- phi_df[which(phi_df[, time_dim] %in% time_laps), ]
        }
    }

    phi_df <- phi_df[, c(kng_dim_1, kng_dim_2, prx)]
    names(phi_df) <- c("src", "dst", "weight")
    grph <- igraph::graph_from_data_frame(phi_df)
    if (directed == FALSE) {
        grph <- igraph::as.undirected(grph)
    }

    w <- igraph::edge.attributes(grph)$weight
    max_span_tree <- igraph::mst(grph,
                                 weights = -w)

    if (!is.null(min_thrashold)) {
        min_thd_subset <- E(grph)$weight > min_thrashold
        min_thd_graph <- subgraph.edges(grph, E(grph)[min_th_subset])
        max_span_tree <- as.undirected(max_span_tree)
        max_span_tree <- union(max_span_tree, min_thd_graph)
    }

    class(max_span_tree) <- c("igraph", "reks_maximum_spanning_tree")

    return(max_span_tree)

}
