#' Simulate RNA-seq Gene Counts
#'
#' Generates a synthetic RNA-seq count matrix from a seed count matrix.
#' If no seed matrix is supplied, the `arab` count matrix bundled with
#' normpatch is used as the seed. The seed matrix is scaled to a common library
#' size and then used as the sampling population by computing row means and
#' variances. Counts are sampled from a negative binomial distribution after
#' applying group-level fold changes.
#'
#' @param n_genes Number of genes to simulate.
#' @param n_replicates Numeric vector giving the number of replicates per group.
#' @param p_DEG Numeric vector giving the proportion of genes assigned as DEGs
#'   for each group.
#' @param fc_mean Numeric vector giving the mean fold change for genes assigned
#'   to each group. When `fc_dist = "fixed"`, this exact value is used.
#' @param fc_dist Distribution used to generate fold changes for genes assigned
#'   as DEGs. One of `"lognormal"`, `"gamma"`, or `"fixed"`. The default is
#'   `"lognormal"`.
#' @param fc_sd Spread parameter for randomized fold changes. For
#'   `fc_dist = "lognormal"`, this is the standard deviation on the natural log
#'   fold change scale. For `fc_dist = "gamma"`, this is the coefficient of
#'   variation.
#' @param group_names Optional group names. If `NULL`, groups are named `"G1"`,
#'   `"G2"`, and so on.
#' @param seed_counts Optional seed count matrix. Rows are genes and columns
#'   are samples. If `NULL`, the package dataset `arab` is used.
#'
#' @return A \linkS4class{SeqCountData} object. Its `meta` slot is a
#'   \linkS4class{SimParams} object.
#'
#' @examples
#' set.seed(1)
#' seed <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#' x <- sim_gene_counts(n_genes = 50, seed_counts = seed)
#' table(def_DEG(x, fc = 4))
#'
#' @export
sim_gene_counts <- function(n_genes = 10000,
                            n_replicates = c(3, 3),
                            p_DEG = c(0.04, 0.01),
                            fc_mean = c(4, 4),
                            fc_dist = c("lognormal", "gamma", "fixed"),
                            fc_sd = 0.3,
                            group_names = NULL,
                            seed_counts = NULL) {
    n_genes <- as.integer(n_genes)
    n_replicates <- as.integer(n_replicates)
    n_groups <- length(n_replicates)
    fc_dist <- match.arg(fc_dist)

    if (length(p_DEG) == 1L) {
        p_DEG <- rep(p_DEG, n_groups)
    }
    if (length(fc_mean) == 1L) {
        fc_mean <- rep(fc_mean, n_groups)
    }
    fc_sd <- rep(as.numeric(fc_sd), length.out = n_groups)

    if (is.null(group_names)) {
        group_names <- paste0("group_", seq_len(n_groups))
    } else {
        group_names <- as.character(group_names)
    }
    seed_counts <- .sim_init_seed_counts(seed_counts)
    population <- .sim_seed_population(seed_counts)

    gene_names <- paste0("gene_", seq_len(n_genes))
    sample_idx <- sample(seq_len(nrow(population)), n_genes, replace = TRUE)
    params <- population[sample_idx, , drop = FALSE]
    rownames(params) <- gene_names

    group_idx <- rep(seq_len(n_groups), times = n_replicates)
    sample_names <- paste0(group_names[group_idx], "_rep", sequence(n_replicates))
    fc <- .sim_fc_matrix(n_genes, p_DEG, fc_mean, fc_dist, fc_sd, group_names)
    rownames(fc) <- gene_names
    sample_fc <- .sim_sample_fc_matrix(fc, group_idx, sample_names)

    data <- apply(
        sample_fc,
        2L,
        function(fc_col, params) {
            stats::rnbinom(
                n = nrow(params),
                mu = fc_col * params$mean,
                size = 1 / params$dispersion
            )
        },
        params = params
    )
    data <- matrix(
        data,
        nrow = n_genes,
        ncol = length(group_idx),
        dimnames = list(gene_names, sample_names)
    )

    exp_design <- data.frame(
        group = factor(group_names[group_idx], levels = group_names),
        row.names = sample_names
    )

    sim_params <- new(
        "SimParams",
        n_genes = as.numeric(n_genes),
        n_groups = as.numeric(n_groups),
        n_replicates = as.numeric(n_replicates),
        p_DEG = as.numeric(p_DEG),
        fc_dist = fc_dist,
        fc_sd = as.numeric(fc_sd),
        fc = fc,
        params_population = population,
        params = params
    )

    new(
        "SeqCountData",
        data = data,
        gene_names = gene_names,
        exp_design = exp_design,
        meta = list(sim_params = sim_params)
    )
}

#' Define Simulated DEGs
#'
#' Uses the fold change matrix stored in `x@meta$sim_params@fc` to define
#' simulated DEG status. A gene is marked `TRUE` when the largest between-group fold
#' change is at least `fc`.
#'
#' @param x A \linkS4class{SeqCountData} object whose `meta` slot is a
#'   \linkS4class{SimParams} object.
#' @param fc Fold change cutoff. Defaults to 2.
#'
#' @return A logical vector with one value per gene.
#'
#' @examples
#' set.seed(1)
#' seed <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#' x <- sim_gene_counts(n_genes = 50, seed_counts = seed)
#' head(def_DEG(x, fc = 4))
#'
#' @export
def_DEG <- function(x, fc = 2) {
    fc_mx <- x@meta$sim_params@fc
    max_fc <- apply(fc_mx, 1L, max)
    min_fc <- apply(fc_mx, 1L, min)
    is_deg <- (max_fc / min_fc) >= fc
    names(is_deg) <- rownames(fc_mx)
    is_deg
}

.sim_init_seed_counts <- function(seed_counts) {
    if (is.null(seed_counts)) {
        data_env <- new.env(parent = emptyenv())
        utils::data("arab", package = "normpatch", envir = data_env)
        seed_counts <- data_env$arab
    }

    seed_counts <- as.matrix(seed_counts)
    storage.mode(seed_counts) <- "numeric"
    seed_counts
}

.sim_seed_population <- function(seed_counts) {
    lib_size <- colSums(seed_counts)
    keep_samples <- is.finite(lib_size) & lib_size > 0
    seed_counts <- seed_counts[, keep_samples, drop = FALSE]
    lib_size <- lib_size[keep_samples]
    seed_counts <- sweep(
        seed_counts,
        2L,
        stats::median(lib_size) / lib_size,
        "*"
    )

    mean_ab <- apply(seed_counts, 1L, mean)
    var_ab <- apply(seed_counts, 1L, stats::var)
    dispersion <- (var_ab - mean_ab) / (mean_ab * mean_ab)
    population <- data.frame(
        mean = mean_ab,
        var = var_ab,
        dispersion = dispersion,
        distribution = "NB"
    )
    keep <- is.finite(population$mean) &
        is.finite(population$var) &
        is.finite(population$dispersion) &
        population$mean > 0 &
        population$dispersion > 0
    population <- population[keep, , drop = FALSE]

    rownames(population) <- paste0("seed_", seq_len(nrow(population)))
    population
}

.sim_fc_matrix <- function(n_genes, p_DEG, fc_mean, fc_dist, fc_sd, group_names) {
    n_groups <- length(group_names)
    fc <- matrix(
        1,
        nrow = n_genes,
        ncol = n_groups,
        dimnames = list(NULL, group_names)
    )
    n_deg <- as.integer(round(n_genes * p_DEG))
    while (sum(n_deg) > n_genes) {
        i <- which.max(n_deg)
        n_deg[i] <- n_deg[i] - 1L
    }

    idx_start <- 1L
    for (i in seq_len(n_groups)) {
        idx_end <- idx_start + n_deg[i] - 1L
        if (idx_end >= idx_start) {
            fc[idx_start:idx_end, i] <- .sim_fc_values(
                n = idx_end - idx_start + 1L,
                mean_fc = fc_mean[i],
                fc_dist = fc_dist,
                fc_sd = fc_sd[i]
            )
            idx_start <- idx_end + 1L
        }
    }
    fc
}

.sim_fc_values <- function(n, mean_fc, fc_dist, fc_sd) {
    if (n < 1L || fc_dist == "fixed" || !is.finite(fc_sd) || fc_sd <= 0) {
        return(rep(mean_fc, n))
    }

    if (fc_dist == "lognormal") {
        return(stats::rlnorm(
            n,
            meanlog = log(mean_fc) - 0.5 * fc_sd * fc_sd,
            sdlog = fc_sd
        ))
    }

    if (fc_dist == "gamma") {
        shape <- 1 / (fc_sd * fc_sd)
        scale <- mean_fc / shape
        return(stats::rgamma(n, shape = shape, scale = scale))
    }

    stop("Unsupported `fc_dist`: ", fc_dist, call. = FALSE)
}

.sim_sample_fc_matrix <- function(fc, group_idx, sample_names) {
    sample_fc <- fc[, group_idx, drop = FALSE]
    colnames(sample_fc) <- sample_names
    sample_fc
}
