skip_if_not_installed("emmeans")

chick_pairs <- function(adjust = "fdr") {
  fit <- lm(weight ~ feed, data = datasets::chickwts)
  emmeans::contrast(emmeans::emmeans(fit, ~ feed), "pairwise",
                    adjust = adjust)
}

pigs_pairs <- function(adjust = "fdr") {
  fit <- lm(log(conc) ~ source * factor(percent), data = emmeans::pigs)
  emmeans::contrast(emmeans::emmeans(fit, ~ source | percent), "pairwise",
                    adjust = adjust)
}

excludes <- function(res, null = 0) {
  res$lower.CL > null | res$upper.CL < null
}

test_that("intervals agree with FDR-adjusted p-values (one family)", {
  res <- confint_fcr(chick_pairs())
  expect_s3_class(res, "summary_emm")
  expect_identical(excludes(res), res$p.value <= 0.05)
})

test_that("the example shows a disagreement in plain emmeans", {
  prs <- chick_pairs()
  bonf <- as.data.frame(confint(prs))
  p_fdr <- as.data.frame(emmeans::test(prs))$p.value
  expect_true(any(p_fdr <= 0.05 & !excludes(bonf)))
})

test_that("p-values equal those of emmeans with adjust = 'fdr'", {
  prs <- chick_pairs()
  res <- confint_fcr(prs)
  expect_equal(res$p.value, as.data.frame(emmeans::test(prs))$p.value)
})

test_that("estimates and standard errors are unchanged", {
  prs <- chick_pairs()
  res <- confint_fcr(prs)
  ref <- as.data.frame(summary(prs))
  expect_equal(res$estimate, ref$estimate)
  expect_equal(res$SE, ref$SE)
})

test_that("intervals are the marginal intervals at the FCR level", {
  prs <- chick_pairs()
  res <- confint_fcr(prs)
  lev <- fcr_level(as.data.frame(emmeans::test(prs, adjust = "none"))$p.value)
  ref <- as.data.frame(confint(prs, level = lev, adjust = "none"))
  expect_equal(res$lower.CL, ref$lower.CL)
  expect_equal(res$upper.CL, ref$upper.CL)
  expect_equal(fcr_info(res)$level, lev)
})

test_that("each 'by' group is a separate family", {
  prs <- pigs_pairs()
  res <- confint_fcr(prs)
  info <- fcr_info(res)
  expect_equal(nrow(info), 4)
  expect_equal(info$m, rep(3, 4))
  expect_equal(res$p.value, as.data.frame(emmeans::test(prs))$p.value)
  expect_identical(excludes(res), res$p.value <= 0.05)

  raw <- as.data.frame(emmeans::test(prs, adjust = "none"))
  for (i in seq_len(nrow(info))) {
    rows <- raw$percent == unique(raw$percent)[i]
    expect_equal(info$level[i], fcr_level(raw$p.value[rows]))
  }
})

test_that("by = NULL pools all contrasts into one family", {
  res <- confint_fcr(pigs_pairs(), by = NULL)
  info <- fcr_info(res)
  expect_equal(nrow(info), 1)
  expect_equal(info$m, 12)
  expect_identical(excludes(res), res$p.value <= 0.05)
})

test_that("ratios on the response scale agree with the p-values", {
  res <- confint_fcr(pigs_pairs(), type = "response")
  expect_true("ratio" %in% names(res))
  expect_identical(excludes(res, null = 1), res$p.value <= 0.05)
})

test_that("odds ratios from a logistic model agree with the p-values", {
  fit <- glm(Pain ~ Treatment + Sex + Age, family = binomial,
             data = emmeans::neuralgia)
  prs <- emmeans::contrast(emmeans::emmeans(fit, ~ Treatment), "pairwise",
                           adjust = "fdr")
  res <- confint_fcr(prs, type = "response")
  expect_true("odds.ratio" %in% names(res))
  lower <- res[[attr(res, "clNames")[1]]]
  upper <- res[[attr(res, "clNames")[2]]]
  expect_identical(lower > 1 | upper < 1, res$p.value <= 0.05)
})

test_that("method = 'BY' is taken from the object or set explicitly", {
  prs_by <- chick_pairs(adjust = "BY")
  res <- confint_fcr(prs_by)
  expect_equal(res$p.value, as.data.frame(emmeans::test(prs_by))$p.value)
  expect_identical(excludes(res), res$p.value <= 0.05)

  res2 <- confint_fcr(chick_pairs(), method = "BY")
  expect_equal(res2$p.value, res$p.value)
  expect_equal(res2$lower.CL, res$lower.CL)
})

test_that("explicit method overrides a non-FDR adjustment", {
  res <- confint_fcr(chick_pairs(adjust = "tukey"), method = "BH")
  expect_equal(res$p.value, confint_fcr(chick_pairs())$p.value)
})

test_that("no rejections gives Bonferroni intervals", {
  fit <- lm(strength ~ machine + diameter, data = emmeans::fiber)
  prs <- emmeans::contrast(emmeans::emmeans(fit, ~ machine), "pairwise",
                           adjust = "fdr")
  res <- confint_fcr(prs)
  bonf <- as.data.frame(confint(prs, adjust = "bonferroni"))
  expect_equal(fcr_info(res)$R, 0)
  expect_equal(res$lower.CL, bonf$lower.CL)
  expect_equal(res$upper.CL, bonf$upper.CL)
  expect_false(any(excludes(res)))
})

test_that("one-sided intervals agree with one-sided p-values", {
  res <- confint_fcr(chick_pairs(), side = ">")
  expect_identical(res$lower.CL > 0, res$p.value <= 0.05)
})

test_that("a different level changes q accordingly", {
  res <- confint_fcr(chick_pairs(), level = 0.90)
  expect_identical(excludes(res), res$p.value <= 0.10)
})

test_that("informative errors", {
  expect_error(confint_fcr(chick_pairs(adjust = "tukey")), "adjust = \"tukey\"")
  expect_error(confint_fcr(chick_pairs(), adjust = "fdr"), "not allowed")
  expect_error(confint_fcr(chick_pairs(), level = 1.5), "strictly between")
  expect_error(confint_fcr(chick_pairs(), method = "holm"), "BH")
  expect_error(confint_fcr(lm(weight ~ feed, data = datasets::chickwts)),
               "No confint_fcr\\(\\) method")
  expect_error(fcr_info(data.frame()), "confint_fcr")
})

test_that("the printed footer describes the adjustment", {
  out <- utils::capture.output(print(confint_fcr(pigs_pairs())))
  expect_true(any(grepl("FCR", out)))
  expect_true(any(grepl("Confidence level used per family", out)))
})

test_that("equivalence and cross-adjust settings cannot slip through", {
  prs <- chick_pairs()
  expect_error(confint_fcr(stats::update(prs, delta = 50)), "equivalence")
  expect_error(confint_fcr(stats::update(pigs_pairs(),
                                           cross.adjust = "bonferroni")),
               "cross.adjust")
  expect_error(confint_fcr(prs, del = 50), "not allowed")
  expect_error(confint_fcr(prs, cross = "bonferroni"), "not allowed")
  expect_error(confint_fcr(prs, adj = "fdr"), "not allowed")
})

test_that("legitimate arguments are not mistaken for forbidden ones", {
  res <- confint_fcr(chick_pairs(), df = 60)
  expect_true(all(res$df == 60))
  expect_identical(excludes(res), res$p.value <= 0.05)
})

test_that("emm_list objects use the contrasts component", {
  fit <- lm(weight ~ feed, data = datasets::chickwts)
  lst <- emmeans::emmeans(fit, pairwise ~ feed, adjust = "fdr")
  res <- confint_fcr(lst)
  expect_equal(res$lower.CL, confint_fcr(chick_pairs())$lower.CL)
  expect_error(confint_fcr(lst, which = "emmeans"), "adjust")
})

test_that("a level stored with update() is used by default", {
  prs <- stats::update(chick_pairs(), level = 0.90)
  res <- confint_fcr(prs)
  expect_identical(excludes(res), res$p.value <= 0.10)
  expect_equal(fcr_info(res)$level,
               fcr_level(as.data.frame(emmeans::test(prs, adjust = "none"))$p.value,
                         q = 0.10))
})

test_that("abbreviated adjust names stored in the object are recognised", {
  expect_equal(fcrci:::resolve_fcr_method(NULL, "fd"), "BH")
  expect_equal(fcrci:::resolve_fcr_method(NULL, "BY"), "BY")
  expect_error(fcrci:::resolve_fcr_method(NULL, "b"), "adjust")
})

test_that("the footer reports the R = 1 convention", {
  fit <- lm(strength ~ machine + diameter, data = emmeans::fiber)
  prs <- emmeans::contrast(emmeans::emmeans(fit, ~ machine), "pairwise",
                           adjust = "BY")
  out <- utils::capture.output(print(confint_fcr(prs)))
  expect_true(any(grepl("R = 1 used", out)))
})
