context("trios - functions for handling id, fid, mid trio information")


## Check output is correct and has not changed
test_that("identify family clusters", {
    id <- 1:11
    fid <- c(NA, NA, 1, 1, NA, 23, 45, 5, 5, 7, NA)
    mid <- c(NA, NA, 2, 2, 65, NA, 46, 6, 6, 6, 0)
    res <- make_family_id(id, fid, mid)
    
    expect_is(res, "integer")
    expect_equal(length(res), length(id))
    expect_equal(res, c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L))
})

