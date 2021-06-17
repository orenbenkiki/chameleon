#' Compute colors for multi-dimensional data.
#'
#' Given a matrix of observation/element rows and variable/measurement columns, compute a color for
#' each row (or group of rows) such that the colors are distinct, and where more-similar colors
#' roughly designate more-similar data rows (or groups of rows).
#'
#' This is intended to provide a "reasonable" set of colors to "arbitrary" data, for use as a
#' convenient default when investigating unknown data sets. It is not meant to replace hand-picked
#' colors tailored for specific data (e.g. using red colors for "bad" rows and green colors for
#' "good" rows).
#'
#' This ensures all colors are distinct by packing the (visible part) of the CIE-L*a*b* color space
#' with the needed number of spheres. To assign the colors to the data, it uses UMAP to reduce the
#' data to 3D. It then uses principal component analysis to represent both the chosen colors (3D
#' sphere centers) and the (3D UMAP) data as point clouds with coordinates in the range 0-1, and
#' finally uses a stable matching algorithm to map these point clouds to each other, thereby
#' assigning a color to each data row. If the data is grouped, then the center of gravity of each
#' group is used to generate a color for each group.
#'
#' By default, "grayscale" colors (with low saturation) are excluded from the result.
#'
#' @param data A matrix whose rows represent elements/observations and columns represent
#'             variables/measurements.
#' @param groups An optional array with an entry per row containing the identifier of the group the
#'               row belongs to.
#' @param umap A boolean specifying whether to run UMAP on the data (by default, `TRUE`). If
#'             `FALSE`, the data matrix must have exactly 3 columns.
#' @param minimal_saturation Exclude colors whose saturation (a^2 + b^2 in CIE-L*a*b* color space)
#'                           is less than this value (by default, 33).
#' @return An array with one entry per row, whose names are the matrix `rownames`, containing the
#'         color of each row. If `groups` was specified, the array will contain one entry per
#'         unique group identifier, whose names are the `as.character` group identifiers, containing
#'         the color of each group.
#'
#' @export
#'
#' @examples
#' chameleon::data_colors(stackloss)
data_colors <- function(data, umap=TRUE, minimal_saturation=33, groups=NULL) {
    if (umap) {
        data <- umap::umap(data, n_components=3)$layout
    }

    if (!is.null(groups)) {
        data <- group_centers(data, groups)
    }

    colors <- pick_n_colors(nrow(data), minimal_saturation)

    data_points <- normalized_data(prcomp(data, retx=TRUE)$x)
    color_points <- normalized_data(prcomp(colors$lab, retx=TRUE)$x)

    distances <- data_distances(data_points, color_points)
    closest_colors <- clue::solve_LSAP(distances)[1:nrow(data)]

    result <- colors$names[closest_colors]
    names(result) <- rownames(data)
    return (result)
}

normalized_data <- function(data) {
    stopifnot(ncol(data) == 3)

    raw_xs <- data[,1]
    raw_ys <- data[,2]
    raw_zs <- data[,3]

    result_xs <- (raw_xs - min(raw_xs)) / (max(raw_xs) - min(raw_xs))
    result_ys <- (raw_ys - min(raw_ys)) / (max(raw_ys) - min(raw_ys))
    result_zs <- (raw_zs - min(raw_zs)) / (max(raw_zs) - min(raw_zs))

    result <- cbind(result_xs, result_ys, result_zs)
    rownames(result) <- rownames(data)
    return (result)
}

group_centers <- function(data, groups) {
    groups <- as.character(groups)
    group_xs <- c()
    group_ys <- c()
    group_zs <- c()
    group_names <- c()
    for (group_name in sort(unique(groups))) {
        group_mask <- groups == group_name
        stopifnot(sum(group_mask) > 0)
        group_x <- mean(na.omit(data[group_mask,1]))
        group_y <- mean(na.omit(data[group_mask,2]))
        group_z <- mean(na.omit(data[group_mask,3]))
        group_xs <- c(group_xs, group_x)
        group_ys <- c(group_ys, group_y)
        group_zs <- c(group_zs, group_z)
        group_names <- c(group_names, group_name)
    }

    result <- cbind(group_xs, group_ys, group_zs)
    rownames(result) <- group_names
    return (result)
}

data_distances <- function(from_data, to_data) {
    distances <- matrix(nrow=nrow(from_data), ncol=nrow(to_data))
    for (from_row in 1:nrow(from_data)) {
        for (to_row in 1:nrow(to_data)) {
            distances[from_row, to_row] <- datum_distance(from_data[from_row,1],
                                                          from_data[from_row,2],
                                                          from_data[from_row,3],
                                                          to_data[to_row,1],
                                                          to_data[to_row,2],
                                                          to_data[to_row,3])
        }
    }

    return (distances)
}

datum_distance <- function(x, y, z, u, v, w) {
    a <- x - u
    b <- y - v
    c <- z - w
    return (sqrt(a * a + b * b + c * c))
}

#' Pick a number of distinct colors which are sufficiently saturated.
#'
#' This ensures all colors are distinct by packing the (visible part) of the CIE-L*a*b* color space
#' with the needed number of spheres and using their centers to generate the colors.
#'
#' @param n The requested (positive) number of colors.
#' @param minimal_saturation Exclude colors whose saturation (`hypot(a, b)` in CIE-L*a*b* color space)
#'                           is less than this value (by default, 33).
#' @return A list with two elements, `names` containing the color names and `lab` containing a
#'         matrix with a row per color and three columns containing the `l`, `a` and `b` coordinates
#'         of each color.
#'
#' @export
#'
#' @examples
#' chameleon::pick_n_colors(8, 33)
pick_n_colors <- function(n, minimal_saturation=33) {
    too_large_step <- 100
    large_step <- 100
    large_step_colors <- pick_step_colors(large_step, minimal_saturation)

    while (is.null(large_step_colors)) {
        too_large_step <- large_step
        large_step <- large_step / 2.0
        large_step_colors <- pick_step_colors(large_step, minimal_saturation)
    }

    if (length(large_step_colors$names) == n) {
        return (large_step_colors)
    }

    if (length(large_step_colors$names) > n) {
        while (TRUE) {
            if (length(large_step_colors$names) == 4) {
                large_step_colors$lab <- head(large_step_colors$lab, n)
                large_step_colors$names <- head(large_step_colors$names, n)
                return (large_step_colors)
            }

            mid_step <- (too_large_step + large_step) / 2
            mid_step_colors <- pick_step_colors(mid_step, minimal_saturation)

            if (is.null(mid_step_colors)) {
                too_large_step <- mid_step
                next
            }

            if (length(mid_step_colors$names) == n) {
                return (mid_step_colors)
            }

            if (length(mid_step_colors$names) < n) {
                large_step <- mid_step
                large_step_colors <- mid_step_colors
            }

            large_step <- mid_step
            large_step_colors <- mid_step_colors
        }
    }

    stopifnot(length(large_step_colors$names) < n)

    small_step <- large_step
    small_step_colors <- large_step_colors
    while (length(small_step_colors$names) < n) {
        small_step <- small_step / 2
        small_step_colors <- pick_step_colors(small_step, minimal_saturation)
        stopifnot(!is.null(small_step_colors))
    }

    stopifnot(length(large_step_colors$names) < n)
    stopifnot(length(small_step_colors$names) >= n)

    while (length(small_step_colors$names) > n) {
        mid_step <- (small_step + large_step) / 2
        mid_step_colors <- pick_step_colors(mid_step, minimal_saturation)
        stopifnot(!is.null(mid_step_colors))
        if (length(mid_step_colors$names) >= n) {
            small_step <- mid_step
            small_step_colors <- mid_step_colors
        } else {
            large_step <- mid_step
            large_step_colors <- mid_step_colors
        }
    }

    stopifnot(length(small_step_colors$names) == n)
    return (small_step_colors)
}

pick_step_colors <- function(step, minimal_saturation) {
    lab <- lab_tetragrid(step)
    srgb <- convertColor(lab, from='Lab', to='sRGB', clip=NA)
    mask <- !is.nan(rowSums(srgb))
    if (sum(mask) < 4) {
        return (NULL)
    }

    lab <- lab[mask,]
    srgb <- srgb[mask,]

    saturation <- sqrt(lab[,2] * lab[,2] + lab[,3] * lab[,3])
    mask <- saturation >= minimal_saturation
    if (sum(mask) < 4) {
        return (NULL)
    }

    lab <- lab[mask,]
    srgb <- srgb[mask,]

    color_names <- as.character(rgb(srgb[,1], srgb[,2], srgb[,3]))

    return (list(lab=lab, names=color_names))
}

lab_tetragrid <- function(step) {
    l_steps <- round(100 / (step * 2 * sqrt(6) / 3))
    ab_steps<- round(250 / step)
    grid <- matrix(nrow=(ab_steps + 1) * (ab_steps + 1) * (l_steps + 1), ncol=3)
    colnames(grid) <- c('l', 'a', 'b')
    i <- 1
    for (la in 0:ab_steps) {
        for (lb in 0:ab_steps) {
            for (li in 0:l_steps) {
                grid[i,3] <- step / 2 - 100 + (2 * la + (lb + li) %% 2) * step
                grid[i,2] <- step / 2 - 150 + (sqrt(3) * (lb + (li %% 2) / 3)) * step
                grid[i,1] <- step / 2 + (li * 2 * sqrt(6) / 3) * step
                i <- i + 1
            }
        }
    }

    return(grid)
}
