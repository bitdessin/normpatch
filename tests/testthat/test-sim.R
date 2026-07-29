make_seed_counts <- function() {
    set.seed(1)
    matrix(
        stats::rnbinom(600, mu = 50, size = 5),
        nrow = 100,
        ncol = 6
    )
}

test_that("sim_gene_counts uses lognormal fold changes by default", {
    x <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.15),
        fc_mean = c(4, 8),
        seed_counts = make_seed_counts()
    )

    fc_group_1 <- x@meta$sim_params@fc[x@meta$sim_params@fc[, 1] != 1, 1]
    fc_group_2 <- x@meta$sim_params@fc[x@meta$sim_params@fc[, 2] != 1, 2]

    expect_identical(x@meta$sim_params@fc_dist, "lognormal")
    expect_equal(x@meta$sim_params@fc_sd, c(0.3, 0.3))
    expect_true(all(fc_group_1 > 0))
    expect_true(all(fc_group_2 > 0))
    expect_gt(length(unique(fc_group_1)), 1)
    expect_gt(length(unique(fc_group_2)), 1)
})

test_that("sim_gene_counts can keep fixed fold changes explicitly", {
    x <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.15),
        fc_mean = c(4, 8),
        fc_dist = "fixed",
        seed_counts = make_seed_counts()
    )

    expect_identical(x@meta$sim_params@fc_dist, "fixed")
    expect_equal(unique(x@meta$sim_params@fc[x@meta$sim_params@fc[, 1] != 1, 1]), 4)
    expect_equal(unique(x@meta$sim_params@fc[x@meta$sim_params@fc[, 2] != 1, 2]), 8)
})

test_that("sim_gene_counts can randomize fold changes", {
    seed <- make_seed_counts()

    set.seed(2)
    x_lognormal <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.25),
        fc_mean = c(4, 4),
        fc_dist = "lognormal",
        fc_sd = 0.5,
        seed_counts = seed
    )
    fc_lognormal <- x_lognormal@meta$sim_params@fc[, 2]
    fc_lognormal <- fc_lognormal[fc_lognormal != 1]

    expect_identical(x_lognormal@meta$sim_params@fc_dist, "lognormal")
    expect_equal(x_lognormal@meta$sim_params@fc_sd, c(0.5, 0.5))
    expect_true(all(fc_lognormal > 0))
    expect_gt(length(unique(fc_lognormal)), 1)

    set.seed(2)
    x_gamma <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.25),
        fc_mean = c(4, 4),
        fc_dist = "gamma",
        fc_sd = 0.5,
        seed_counts = seed
    )
    fc_gamma <- x_gamma@meta$sim_params@fc[, 2]
    fc_gamma <- fc_gamma[fc_gamma != 1]

    expect_identical(x_gamma@meta$sim_params@fc_dist, "gamma")
    expect_true(all(fc_gamma > 0))
    expect_gt(length(unique(fc_gamma)), 1)
})
