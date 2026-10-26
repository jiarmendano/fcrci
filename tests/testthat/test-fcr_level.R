test_that("fcr_level matches the closed-form BH level", {
  p <- c(0.0004, 0.003, 0.012, 0.04, 0.21, 0.67)
  n_rej <- sum(p.adjust(p, "BH") <= 0.05)
  expect_equal(n_rej, 3)
  expect_equal(fcr_level(p, q = 0.05), 1 - n_rej * 0.05 / 6)
})

test_that("fcr_level matches the closed-form BY level", {
  p <- c(0.0004, 0.003, 0.012, 0.04, 0.21, 0.67)
  c_m <- sum(1 / 1:6)
  n_rej <- sum(p.adjust(p, "BY") <= 0.05)
  expect_equal(fcr_level(p, q = 0.05, method = "BY"),
               1 - n_rej * 0.05 / (6 * c_m))
})

test_that("'fdr' is accepted as a synonym of 'BH'", {
  p <- c(0.0004, 0.003, 0.012, 0.04, 0.21, 0.67)
  expect_equal(fcr_level(p, method = "fdr"), fcr_level(p, method = "BH"))
})

test_that("no rejections gives the Bonferroni level", {
  p <- c(0.2, 0.5, 0.9)
  expect_equal(fcr_level(p, q = 0.05), 1 - 0.05 / 3)
  expect_equal(fcr_level(p, q = 0.05, method = "BY"),
               1 - 0.05 / (3 * sum(1 / 1:3)))
})

test_that("all rejections gives the unadjusted level", {
  expect_equal(fcr_level(c(0.001, 0.002, 0.01), q = 0.05), 0.95)
})

test_that("missing p-values are ignored", {
  expect_equal(fcr_level(c(0.001, NA, 0.5, 0.2), q = 0.05),
               fcr_level(c(0.001, 0.5, 0.2), q = 0.05))
})

test_that("intervals at the FCR level agree with adjusted p-values", {
  set.seed(1)
  for (i in 1:200) {
    m <- sample(2:40, 1)
    p <- c(runif(m) ^ sample(c(1, 3, 6), 1))
    for (meth in c("BH", "BY")) {
      lev <- fcr_level(p, q = 0.05, method = meth)
      expect_identical(p <= 1 - lev, p.adjust(p, meth) <= 0.05)
    }
  }
})

test_that("fcr_level validates its input", {
  expect_error(fcr_level("a"), "numeric")
  expect_error(fcr_level(c(0.1, 0.2), q = 0), "strictly between")
  expect_error(fcr_level(c(0.1, 0.2), q = c(0.05, 0.1)), "strictly between")
  expect_error(fcr_level(c(NA_real_, NA_real_)), "no non-missing")
  expect_error(fcr_level(c(0.1, 1.2)), "\\[0, 1\\]")
  expect_error(fcr_level(c(0.1, 0.2), method = "holm"))
})
