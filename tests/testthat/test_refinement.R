

test_that("refinement hyper category test works fine",{
    expect_true(all.equal(hyper(c(3,4,6), c(12,12,12), c(0,2,12), c(19,19,19)), c(0.0489433, 0.1366560, 0.863453), tolerance=1.5e-6))
    expect_true(all.equal(hyper(c(3,4,6), c(12,12,12), c(0,2,12), c(19,19,19), under=T), c(1, 0.978307, 0.362281), tolerance=1.5e-6))

})

