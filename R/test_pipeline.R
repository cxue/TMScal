library(testthat)
library(TMScal)

test_that("generate_pathway_names returns 192 pathways", {
    pathways <- generate_pathway_names()
    expect_equal(length(pathways), 192)
})

test_that("get_pathway works with 11-nt context", {
    result <- get_pathway("C", "T", "TGTTCGAGCTC")
    expect_equal(result, "CGA_C_T")
})

test_that("get_pathway works with 3-nt context", {
    result <- get_pathway("C", "T", "CGA")
    expect_equal(result, "CGA_C_T")
})

test_that("stratified_split maintains event rate", {
    clin <- data.frame(
        dfs_status = c(rep(0, 70), rep(1, 30)),
        dfs_time = runif(100, 1, 100)
    )
    split <- stratified_split(clin, 0.6, 42)
    train_rate <- mean(clin$dfs_status[split$train])
    test_rate <- mean(clin$dfs_status[split$test])
    expect_equal(train_rate, test_rate, tolerance = 0.05)
})