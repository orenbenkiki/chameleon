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
#' This ensures all colors are distinct by packing the (visible part) of the CIELAB color space
#' with the needed number of spheres. To assign the colors to the data, it uses UMAP to reduce the
#' data to 3D. It then uses principal component analysis to represent both the chosen colors (3D
#' sphere centers) and the (3D UMAP) data as point clouds with coordinates in the range 0-1, and
#' finally uses a stable matching algorithm to map these point clouds to each other, thereby
#' assigning a color to each data row. If the data is grouped, then the center of gravity of each
#' group is used to generate a color for each group.
#'
#' @param data A matrix whose rows represent elements/observations and columns represent
#'             variables/measurements.
#' @param groups An optional array with an entry per row containing the identifier of the group the
#'               row belongs to.
#' @param run_umap A boolean specifying whether to run UMAP on the data to convert it to 3D (by default,
#'                 \code{TRUE}). If \code{FALSE}, the data matrix must have exactly 3 columns and
#'                 will be used as-is.
#' @param minimal_saturation Exclude colors whose saturation (\code{hypot(a, b)} in CIELAB color
#'                           space) is less than this value (by default, 33).
#' @param minimal_lightness Exclude colors whose lightnes (\code{l} in CIELAB color space) is less
#'                          than this value (by default, 20).
#' @param maximal_lightness Exclude colors whose lightnes (\code{l} in CIELAB color space) is more
#'                          than this value (by default, 80).
#' @return An array with one entry per row, whose names are the matrix \code{rownames}, containing the
#'         color of each row. If \code{groups} was specified, the array will contain one entry per
#'         unique group identifier, whose names are the \code{as.character} group identifiers,
#'         containing the color of each group.
#'
#' @export
#'
#' @examples
#' chameleon::data_colors(stackloss)
data_colors <- function(data, run_umap=TRUE, groups=NULL,
                        minimal_saturation=33, minimal_lightness=20, maximal_lightness=80) {
    if (run_umap) {
        data <- umap::umap(data, n_components=3)$layout
    }

    if (!is.null(groups)) {
        data <- group_centers(data, groups)
    }

    colors <- distinct_colors(nrow(data), minimal_saturation, minimal_lightness, maximal_lightness)

    data_points <- normalized_data(stats::prcomp(data, retx=TRUE)$x)
    color_points <- normalized_data(stats::prcomp(colors$lab, retx=TRUE)$x)

    distances <- data_distances(data_points, color_points)
    closest_colors <- clue::solve_LSAP(distances)[1:nrow(data)]

    result <- colors$name[closest_colors]
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
        group_x <- mean(stats::na.omit(data[group_mask,1]))
        group_y <- mean(stats::na.omit(data[group_mask,2]))
        group_z <- mean(stats::na.omit(data[group_mask,3]))
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

#' Pick a number of distinct colors.
#'
#' This ensures all colors are distinct by packing the (visible part) of the CIELAB color space
#' with the needed number of spheres, and using their centers to generate the colors.
#'
#' @param n The requested (positive) number of colors.
#' @param minimal_saturation Exclude colors whose saturation (\code{hypot(a, b)} in CIELAB color
#'                           space) is less than this value (by default, 33).
#' @param minimal_lightness Exclude colors whose lightnes (\code{l} in CIELAB color space) is less
#'                          than this value (by default, 20).
#' @param maximal_lightness Exclude colors whose lightnes (\code{l} in CIELAB color space) is more
#'                          than this value (by default, 80).
#' @return A list with two elements, \code{name} containing the color names and \code{lab}
#'         containing a matrix with a row per color and three columns containing the \code{l},
#'         \code{a} and \code{b} coordinates of each color.
#'
#' @export
#'
#' @examples
#' chameleon::distinct_colors(8)
distinct_colors <- function(n, minimal_saturation=33, minimal_lightness=20, maximal_lightness=80) {
    stopifnot(n > 0)

    stopifnot(minimal_saturation > 0)
    stopifnot(minimal_saturation < 150)

    stopifnot(minimal_lightness >= 0)
    stopifnot(minimal_lightness <= maximal_lightness)
    stopifnot(maximal_lightness <= 100)

    too_large_step <- 100
    large_step <- 100
    large_step_colors <- pick_step_colors(large_step,
                                          minimal_saturation,
                                          minimal_lightness,
                                          maximal_lightness)

    while (is.null(large_step_colors)) {
        too_large_step <- large_step
        large_step <- large_step / 2.0
        large_step_colors <- pick_step_colors(large_step,
                                              minimal_saturation,
                                              minimal_lightness,
                                              maximal_lightness)
    }

    if (length(large_step_colors$name) == n) {
        return (large_step_colors)
    }

    if (length(large_step_colors$name) > n) {
        while (TRUE) {
            if (length(large_step_colors$name) == 4) {
                large_step_colors$lab <- utils::head(large_step_colors$lab, n)
                large_step_colors$name <- utils::head(large_step_colors$name, n)
                return (large_step_colors)
            }

            mid_step <- (too_large_step + large_step) / 2
            mid_step_colors <- pick_step_colors(mid_step,
                                                minimal_saturation,
                                                minimal_saturation,
                                                maximal_lightness)

            if (is.null(mid_step_colors)) {
                too_large_step <- mid_step
                next
            }

            if (length(mid_step_colors$name) == n) {
                return (mid_step_colors)
            }

            if (length(mid_step_colors$name) < n) {
                large_step <- mid_step
                large_step_colors <- mid_step_colors
            }

            large_step <- mid_step
            large_step_colors <- mid_step_colors
        }
    }

    stopifnot(length(large_step_colors$name) < n)

    small_step <- large_step
    small_step_colors <- large_step_colors
    while (length(small_step_colors$name) < n) {
        small_step <- small_step / 2
        small_step_colors <- pick_step_colors(small_step,
                                              minimal_saturation,
                                              minimal_lightness,
                                              maximal_lightness)
        stopifnot(!is.null(small_step_colors))
    }

    stopifnot(length(large_step_colors$name) < n)
    stopifnot(length(small_step_colors$name) >= n)

    while (length(small_step_colors$name) > n) {
        if (large_step - small_step < 1e-6) {
            small_step_colors$name <- utils::head(small_step_colors$name, n)
            small_step_colors$lab <- utils::head(small_step_colors$lab, n)
            break
        }

        mid_step <- (small_step + large_step) / 2
        mid_step_colors <- pick_step_colors(mid_step,
                                            minimal_saturation,
                                            minimal_lightness,
                                            maximal_lightness)
        stopifnot(!is.null(mid_step_colors))
        if (length(mid_step_colors$name) >= n) {
            small_step <- mid_step
            small_step_colors <- mid_step_colors
        } else {
            large_step <- mid_step
            large_step_colors <- mid_step_colors
        }
    }

    stopifnot(length(small_step_colors$name) == n)
    return (small_step_colors)
}

pick_step_colors <- function(step, minimal_saturation, minimal_lightness, maximal_lightness) {
    lab <- lab_tetragrid(step)
    srgb <- grDevices::convertColor(lab, from='Lab', to='sRGB', clip=NA)
    mask <- !is.nan(rowSums(srgb))
    if (sum(mask) < 4) {
        return (NULL)
    }

    lab <- lab[mask,]
    srgb <- srgb[mask,]

    saturation <- sqrt(lab[,2] * lab[,2] + lab[,3] * lab[,3])
    mask <- (saturation >= minimal_saturation) & (lab[,1] >= minimal_lightness) & (lab[,1] <= maximal_lightness)
    if (sum(mask) < 4) {
        return (NULL)
    }

    lab <- lab[mask,]
    srgb <- srgb[mask,]

    color_names <- as.character(grDevices::rgb(srgb[,1], srgb[,2], srgb[,3]))

    return (list(lab=lab, name=color_names))
}

lab_tetragrid <- function(step) {
    l_steps <- round(100 / (step * 2 * sqrt(6) / 3))
    a_steps<- round(250 / (step * 2))
    b_steps<- round(250 / (step * sqrt(3)))
    grid <- matrix(nrow=(l_steps + 1) * (a_steps + 1) * (b_steps + 1), ncol=3)
    colnames(grid) <- c('l', 'a', 'b')
    i <- 1
    for (li in 0:l_steps) {
        for (la in 0:a_steps) {
            for (lb in 0:b_steps) {
                grid[i,1] <- step / 2 + (li * 2 * sqrt(6) / 3) * step
                grid[i,2] <- step / 2 - 100 + (2 * la + (lb + li) %% 2) * step
                grid[i,3] <- step / 2 - 150 + (sqrt(3) * (lb + (li %% 2) / 3)) * step
                i <- i + 1
            }
        }
    }

    return(grid)
}

#' Setup a color scale of distinct discrete colors in ggplot2.
#'
#' This is a thin wrapper to \code{ggplot2::discrete_scale('colour', 'chameleon', ...)}, which uses
#' the colors chosen by invoking \code{distinct_colors}. The order of the colors is arbitrary. If
#' the data has some structure the colors should reflect, use one of the many palettes available in
#' R, or using \code{data_colors} for automatically matching the colors to the structure of
#' multi-dimensional data.
#'
#' @param ... Additional parameters for \code{discrete_scale}.
#'
#' @inheritParams distinct_colors
#'
#' @examples
#' library(ggplot2)
#' data(pbmc)
#' frame <- as.data.frame(pbmc$umap)
#' frame$type <- pbmc$types
#' ggplot(frame, aes(x=xs, y=ys, color=type)) +
#'     geom_point(size=0.75) +
#'     scale_color_chameleon() +
#'     theme(legend.text=element_text(size=12), legend.key.height=unit(14, 'pt'))
#' @export
scale_color_chameleon <- function(minimal_saturation = 33, minimal_lightness = 20, maximal_lightness = 80, ...) {
    color_func <- function(n) distinct_colors(n, minimal_saturation = minimal_saturation,
                                              minimal_lightness = minimal_lightness,
                                              maximal_lightness = maximal_lightness)$name
    ggplot2::discrete_scale('colour', 'chameleon', color_func, ...)
}

#' Setup a fill scale of distinct discrete colors in ggplot2.
#'
#' This is a thin wrapper to \code{ggplot2::discrete_scale('fill', 'chameleon', ...)}, which uses
#' the colors chosen by invoking \code{distinct_colors}. The order of the colors is arbitrary. If
#' the data has some structure the colors should reflect, use one of the many palettes available in
#' R, or using \code{data_colors} for automatically matching the colors to the structure of
#' multi-dimensional data.
#'
#' @param ... Additional parameters for \code{discrete_scale}.
#'
#' @inheritParams distinct_colors
#'
#' @examples
#' library(ggplot2)
#' data(pbmc)
#' frame <- as.data.frame(pbmc$umap)
#' frame$type <- pbmc$types
#' ggplot(frame, aes(x=xs, y=ys, fill=type)) +
#'     geom_point(size=0.75, shape=21, color="black", stroke=0.1) +
#'     scale_fill_chameleon() +
#'     theme(legend.text=element_text(size=12), legend.key.height=unit(14, 'pt'))
#' @export
scale_fill_chameleon <- function(minimal_saturation = 33, minimal_lightness = 20, maximal_lightness = 80, ...) {
    color_func <- function(n) distinct_colors(n, minimal_saturation = minimal_saturation,
                                              minimal_lightness = minimal_lightness,
                                              maximal_lightness = maximal_lightness)$name
    ggplot2::discrete_scale('fill', 'chameleon', color_func, ...)
}
