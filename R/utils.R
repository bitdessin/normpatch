#' Get Normalized Counts
#'
#' Returns normalized counts using effective library sizes.
#'
#' @param x A `SeqCountData` object.
#' @return A normalized count matrix.
#'
#' @examples
#' counts <- matrix(1:12, nrow = 3)
#' colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#' rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#' exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#' x <- newSeqCountData(
#'     counts,
#'     exp_design = exp_design,
#'     norm_factors = rep(1, ncol(counts))
#' )
#' norm_counts(x)
#'
#' @export
norm_counts <- function(x) {
    if (is.null(x@meta$nf)) {
        stop("Normalization factors are not set. Run `calc_nf()` first.", call. = FALSE)
    }
    sweep(x@data, 2, .calc_sizefactors(x@data, x@meta$nf), "/")
}

#' Plot an MA Plot
#'
#' Creates an MA plot from normalized counts for two groups in `x@exp_design$group`.
#' Group averages are computed before log transformation for the A values.
#'
#' @param x A `SeqCountData` object.
#' @param comparison Optional vector of two group names to compare. If
#'   `NULL`, the first two groups in `x@exp_design$group` are used.
#' @param esp Value added before log transformation.
#' @param col Optional point colors. Can be a single color or one color per gene.
#'   When provided, points are drawn from the most frequent color category to
#'   the least frequent category.
#' @return A `ggplot` object.
#' @importFrom ggplot2 aes geom_point ggplot labs scale_color_identity theme_bw
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(
#'         counts,
#'         exp_design = exp_design,
#'         norm_factors = rep(1, ncol(counts))
#'     )
#'     plot_ma(x)
#' }
#'
#' @export
plot_ma <- function(x, comparison = NULL, esp = 1, col = NULL) {
    .require_namespace("ggplot2")

    group_vec <- factor(x@exp_design$group)
    groups <- levels(group_vec)
    if (length(groups) < 2L) {
        stop("At least two groups are required for an MA plot.", call. = FALSE)
    }
    if (is.null(comparison)) {
        comparison <- groups[seq_len(2L)]
    } else {
        if (length(comparison) != 2L) {
            stop("`comparison` must be NULL or a length-2 vector.", call. = FALSE)
        }
        comparison <- as.character(comparison)
        missing_groups <- setdiff(comparison, groups)
        if (length(missing_groups) > 0L) {
            stop(
                "`comparison` contains groups not present in `x@exp_design$group`: ",
                paste(missing_groups, collapse = ", "),
                call. = FALSE
            )
        }
    }
    data <- norm_counts(x)

    y1 <- rowMeans(data[, group_vec == comparison[1], drop = FALSE])
    y2 <- rowMeans(data[, group_vec == comparison[2], drop = FALSE])
    log_y1 <- log2(y1 + esp)
    log_y2 <- log2(y2 + esp)

    data_df <- data.frame(
        A = unname(log2(((y1 + y2) / 2) + esp)),
        M = unname(log_y2 - log_y1),
        id = as.character(.get_gene_names(x)),
        check.names = FALSE
    )

    if (!is.null(col)) {
        data_df$col <- col
        col_counts <- sort(table(data_df$col), decreasing = TRUE)
        data_df <- data_df[
            order(match(data_df$col, names(col_counts))),
            ,
            drop = FALSE
        ]
    }

    p <- ggplot(data_df, aes(x = .data[["A"]], y = .data[["M"]], text = .data[["id"]]))
    if (is.null(col)) {
        p <- p + geom_point(alpha = 0.4, size = 1.2)
    } else {
        p <- p +
            geom_point(
                aes(color = .data[["col"]]),
                alpha = 0.4,
                size = 1.2
            ) +
            scale_color_identity()
    }

    p +
        labs(
            x = paste0('log2((', comparison[1], ' + ', comparison[2], ') / 2)'),
            y = paste0('log2(', comparison[2], ') - log2(', comparison[1], ')'),
            title = paste(comparison[2], "vs", comparison[1])
        ) +
        theme_bw()
}

#' Plot a PCA Projection
#'
#' Creates a PCA plot from log-normalized counts.
#'
#' @param x A `SeqCountData` object.
#' @param pc_x Principal component to show on the horizontal axis.
#' @param pc_y Principal component to show on the vertical axis.
#' @param esp Value added before log transformation.
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(
#'         counts,
#'         exp_design = exp_design,
#'         norm_factors = rep(1, ncol(counts))
#'     )
#'     plot_pca(x)
#' }
#'
#' @export
plot_pca <- function(x, pc_x = 1, pc_y = 2, esp = 1) {
    .require_namespace("ggplot2")

    data <- log2(norm_counts(x) + esp)
    max_pc <- min(nrow(data), ncol(data))
    if (max_pc < 2L) {
        stop("At least two genes and two samples are required for a PCA plot.", call. = FALSE)
    }
    if (pc_x < 1L || pc_y < 1L || pc_x > max_pc || pc_y > max_pc) {
        stop("`pc_x` and `pc_y` must refer to available principal components.", call. = FALSE)
    }

    pca <- stats::prcomp(t(data), center = TRUE, scale. = FALSE)
    variance <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

    data_df <- data.frame(
        sample_id = colnames(data),
        x = pca$x[, pc_x],
        y = pca$x[, pc_y],
        group = x@exp_design$group,
        check.names = FALSE
    )

    ggplot(
        data_df,
        aes(
            x = .data[["x"]],
            y = .data[["y"]],
            color = .data[["group"]]
        )
    ) +
        geom_point(size = 2.4) +
        labs(
            x = paste0("PC", pc_x, " (", round(variance[pc_x], 1), "%)"),
            y = paste0("PC", pc_y, " (", round(variance[pc_y], 1), "%)"),
            color = "group"
        ) +
        theme_bw()
}


.require_namespace <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Package `", package, "` is required for this function.", call. = FALSE)
    }
    invisible(TRUE)
}

.calc_sizefactors <- function(x, nf) {
    effective_libsizes <- nf * colSums(x)
    effective_libsizes / mean(effective_libsizes)
}

.get_gene_names <- function(x) {
    if (!is.null(x@gene_names)) {
        return(x@gene_names)
    }
    if (!is.null(rownames(x@data))) {
        return(rownames(x@data))
    }
    paste0("gene_", seq_len(nrow(x@data)))
}
