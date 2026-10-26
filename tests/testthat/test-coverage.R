test_that("FCR is controlled for independent estimates (simulation)", {
  skip_on_cran()
  set.seed(2005)
  m <- 40
  q <- 0.05
  theta <- c(rep(3, 10), rep(0, m - 10))
  n_rep <- 3000
  fcp_fcr <- fcp_naive <- numeric(n_rep)
  for (b in seq_len(n_rep)) {
    est <- rnorm(m, theta)
    p <- 2 * pnorm(-abs(est))
    selected <- p.adjust(p, "BH") <= q
    if (!any(selected)) next
    half_fcr <- qnorm(1 - (1 - fcr_level(p, q)) / 2)
    half_naive <- qnorm(1 - q / 2)
    miss_fcr <- abs(est - theta) > half_fcr
    miss_naive <- abs(est - theta) > half_naive
    fcp_fcr[b] <- mean(miss_fcr[selected])
    fcp_naive[b] <- mean(miss_naive[selected])
  }
  mc_se <- stats::sd(fcp_fcr) / sqrt(n_rep)
  expect_lte(mean(fcp_fcr), q + 3 * mc_se)
  expect_gt(mean(fcp_naive), mean(fcp_fcr))
})
