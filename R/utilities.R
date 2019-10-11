# @importFrom methods is
#' @importFrom methods as

rta <- function (data, binary = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop(paste0("Package \"Matrix\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }
    spm <- any(methods::is(data) == "sparseMatrix")
    num <- Matrix::t(data / Matrix::rowSums(data))
    den <- Matrix::colSums(data) / sum(data)
    RA <- Matrix::t(num / den)
    RA[is.na(RA)] <- 0
    if (isTRUE(binary)) {
        if (isTRUE(spm)) {
            RA <- as(as(RA >= 1, "lgCMatrix"), "dgCMatrix")
        } else {
            RA <- ifelse(RA >= 1, 1, 0)
            RA <- Matrix::Matrix(RA)
        }
    }
    return(RA)
}

#' @importFrom stats setNames reshape

wide_to_long <- function(df, n_dims, col_names, measure_name) {
    # Create long table
    if (n_dims == 3) {
        if (is.matrix(df)) {
            rnms <- colnames(df)
            df <- cbind.data.frame(rownames(df), as.data.frame(df))
            names(df) <- c(col_names[1], rnms)
        } else {
            rnms <- names(df)
            df <- Map(function(x, y) setNames(x, c(col_names[1], y)),
                      df, rnms)
            df <- Reduce(function(x1, x2)
                merge(x1, x2, all = TRUE, by = col_names[1]), df)
        }
        df <- reshape(df, idvar = col_names[1], varying = list(2:ncol(df)),
                      v.names = measure_name, times = rnms,
                      direction = "long")
    } else {
        df <- cbind.data.frame(names(df), as.vector(df),
                               stringsAsFactors = FALSE)
    }
    # Assign the name of the measure
    colnames(df) <- c(col_names, measure_name)
    # Transform the temporal dimension in a numeric variable
    if (is.factor(df[, 2])) {
        df[, 2] <- as.numeric(levels(df[, 2]))[df[, 2]]
    } else {
        df[, 2] <- as.numeric(df[, 2])
    }
    # Order the data.frame by location name (and time)
    if (n_dims == 3) {
        ord <- order(as.character(df[, 2]),
                     as.character(df[, 1]))
    } else {
        ord <- order(as.character(df[, 1]))
    }

    df <- df[ord, ]
    rownames(df) <- NULL
    return(df)
}

data_info <- function(occt) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop(paste0("Package \"Matrix\" needed for this function to work. ",
                    "Please install it."), call. = FALSE)
    }
    # if (!"table" %in% class(occt)) {
    #     stop("Please, provide an object of class \"table\" as input.")
    # }
    dims <- dim(occt)
    n_dims <- length(dims)
    dim_nms <- names(attr(occt, "dimnames"))
    if (!n_dims %in% 2:3) {
        stop(paste0("You can apply the function to ",
                    "a 2d or 3d object, ",
                    "nothing more, nothing less."))
    }
    if (n_dims == 3) {
        nd1 <- c(1, 3)
        nd2 <- c(2, 3)
    } else {
        nd1 <- 1
        nd2 <- 2
    }
    return(list(n_dims = n_dims,
                dim_nms = dim_nms,
                nd = list(nd1, nd2)))
}

remove_zeros <- function(occt) {
    if (any(Matrix::rowSums(occt) == 0)) {
        occt <- occt[-Matrix::which(Matrix::rowSums(occt) == 0), ]
    }
    if (any(Matrix::colSums(occt) == 0)) {
        occt <- occt[, -Matrix::which(Matrix::colSums(occt) == 0)]
    }
    return(occt)
}
