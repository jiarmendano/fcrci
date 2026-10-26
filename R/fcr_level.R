#' Confidence level for FCR-adjusted confidence intervals
#'
#' Computes the confidence level that makes marginal confidence intervals
#' consistent with Benjamini-Hochberg (BH) or Benjamini-Yekutieli (BY) adjusted
#' p-values, following the false coverage-statement rate (FCR) adjustment of
#' Benjamini and Yekutieli (2005).
#'
#' Let \eqn{m} be the number of non-missing p-values in the family and \eqn{R}
#' the number of hypotheses rejected at FDR level \eqn{q}. The FCR-adjusted
#' confidence level is
#' \deqn{1 - \frac{R \, q}{m \, c_m},}{1 - R * q / (m * c_m),}
#' where \eqn{c_m = 1} for `method = "BH"` and
#' \eqn{c_m = \sum_{i=1}^{m} 1/i}{c_m = sum(1 / (1:m))} for `method = "BY"`.
#'
#' Because the step-up procedure rejects exactly the hypotheses with
#' \eqn{p_i \le R q / (m c_m)}{p_i <= R * q / (m * c_m)}, a marginal interval
#' computed at this level excludes the null value if and only if the
#' corresponding adjusted p-value is at most \eqn{q}, provided that the
#' interval and the unadjusted p-value come from the same test statistic and
#' degrees of freedom.
#'
#' When no hypothesis is rejected (\eqn{R = 0}), the level is computed with
#' \eqn{R = 1}, that is \eqn{1 - q / (m c_m)}: the Bonferroni level for
#' `"BH"`, and the Bonferroni level at \eqn{q / c_m} for `"BY"`. Benjamini and
#' Yekutieli (2005) construct no intervals when nothing is selected; this
#' convention of the package keeps intervals and adjusted p-values consistent,
#' because no interval then excludes the null value.
#'
#' @param p Numeric vector of unadjusted p-values forming one family of tests.
#'   Missing values are ignored.
#' @param q FCR (and FDR) level, a single number strictly between 0 and 1.
#'   Usually `1 - level`, for example 0.05 for 95% intervals.
#' @param method Either `"BH"` (Benjamini-Hochberg selection; `"fdr"`, as in
#'   [stats::p.adjust()] and 'emmeans', is accepted as a synonym) or `"BY"`
#'   (Benjamini-Yekutieli selection, valid under arbitrary dependence).
#'
#' @return A single number: the confidence level to use for the marginal
#'   intervals of the family.
#'
#' @references
#' Benjamini, Y. and Yekutieli, D. (2005). False discovery rate-adjusted
#' multiple confidence intervals for selected parameters. *Journal of the
#' American Statistical Association*, 100(469), 71-81.
#' \doi{10.1198/016214504000001907}
#'
#' Jung, K., Friede, T. and Beissbarth, T. (2011). Reporting FDR analogous
#' confidence intervals for the log fold change of differentially expressed
#' genes. *BMC Bioinformatics*, 12, 288. \doi{10.1186/1471-2105-12-288}
#'
#' @seealso [confint_fcr()] to obtain the intervals for 'emmeans' contrasts.
#'
#' @examples
#' p <- c(0.0004, 0.003, 0.012, 0.04, 0.21, 0.67)
#' lev <- fcr_level(p, q = 0.05)
#' lev
#'
#' # The level is consistent with the BH-adjusted p-values:
#' # an interval at level `lev` excludes the null exactly when p <= 1 - lev.
#' p.adjust(p, method = "BH") <= 0.05
#' p <= 1 - lev
#'
#' # Valid under arbitrary dependence, at the cost of wider intervals
#' fcr_level(p, q = 0.05, method = "BY")
#' @export
fcr_level <- function(p, q = 0.05, method = c("BH", "BY")) {
  if (is.character(method) && length(method) == 1L &&
      toupper(method) == "FDR") {
    method <- "BH"
  }
  method <- match.arg(method)
  fcr_family(p, q, method)$level
}

# Internal workhorse: returns level, m and R for one family.
fcr_family <- function(p, q, method) {
  if (!is.numeric(p)) {
    stop("'p' must be a numeric vector of p-values.", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0 || q >= 1) {
    stop("'q' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  p <- p[!is.na(p)]
  if (length(p) == 0L) {
    stop("'p' contains no non-missing p-values.", call. = FALSE)
  }
  if (any(p < 0 | p > 1)) {
    stop("All p-values must lie in [0, 1].", call. = FALSE)
  }

  m <- length(p)
  n_rejected <- sum(stats::p.adjust(p, method = method) <= q)
  c_m <- if (method == "BY") sum(1 / seq_len(m)) else 1
  level <- 1 - max(n_rejected, 1L) * q / (m * c_m)

  list(level = level, m = m, R = n_rejected)
}
