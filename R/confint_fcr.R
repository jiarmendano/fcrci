#' FDR-consistent (FCR-adjusted) confidence intervals
#'
#' Computes confidence intervals that agree with false discovery rate (FDR)
#' adjusted p-values, using the false coverage-statement rate (FCR) adjustment
#' of Benjamini and Yekutieli (2005). Within each family of tests, an interval
#' excludes the null value if and only if the corresponding FDR-adjusted
#' p-value is at most `1 - level`.
#'
#' `confint_fcr()` is a generic function. The current release provides methods
#' for `emmGrid` and `emm_list` objects created by the 'emmeans' package,
#' typically the result of [emmeans::contrast()] or `pairs()`.
#'
#' @param object An object for which FCR-adjusted intervals are wanted.
#' @param ... Further arguments passed to methods.
#'
#' @return See the documentation of the specific method, for example
#'   [confint_fcr.emmGrid()].
#'
#' @references
#' Benjamini, Y. and Yekutieli, D. (2005). False discovery rate-adjusted
#' multiple confidence intervals for selected parameters. *Journal of the
#' American Statistical Association*, 100(469), 71-81.
#' \doi{10.1198/016214504000001907}
#'
#' @seealso [confint_fcr.emmGrid()], [fcr_level()], [fcr_info()].
#' @export
confint_fcr <- function(object, ...) {
  UseMethod("confint_fcr")
}

#' @rdname confint_fcr
#' @export
confint_fcr.default <- function(object, ...) {
  stop("No confint_fcr() method for objects of class '",
       paste(class(object), collapse = "/"), "'. ",
       "Supported classes: 'emmGrid' and 'emm_list' (package 'emmeans').",
       call. = FALSE)
}

#' FCR-adjusted confidence intervals for 'emmeans' contrasts
#'
#' Returns the usual 'emmeans' summary table (estimates, confidence intervals,
#' test statistics and p-values) in which the p-values are FDR-adjusted and the
#' confidence intervals are FCR-adjusted, so that both lead to the same
#' conclusions.
#'
#' @section Why this is needed:
#' In 'emmeans', the `adjust` argument accepts FDR methods such as `"fdr"`
#' (`"BH"`) or `"BY"` for p-values, but confidence intervals cannot be
#' FDR-adjusted, so `confint()` switches to a Bonferroni adjustment (noted only
#' in the footer of the output). As a result, a contrast can be significant
#' according to its FDR-adjusted p-value while its (Bonferroni-adjusted)
#' confidence interval includes the null value. `confint_fcr()` replaces those
#' intervals with the FCR-adjusted intervals of Benjamini and Yekutieli (2005),
#' which are the interval counterpart of the Benjamini-Hochberg procedure.
#'
#' @section Method:
#' For each family of \eqn{m} tests, the FDR procedure rejects \eqn{R}
#' hypotheses at level \eqn{q = 1 - \mathrm{level}}{q = 1 - level}. Marginal
#' intervals are then computed at confidence level \eqn{1 - R q / m} (or
#' \eqn{1 - R q / (m \sum_{i=1}^m 1/i)}{1 - R q / (m sum(1/i))} for
#' `method = "BY"`); see [fcr_level()]. When \eqn{R = 0}, the level is
#' computed with \eqn{R = 1}. Benjamini and Yekutieli (2005) construct no
#' intervals when nothing is selected; using \eqn{R = 1} is a convention of
#' this package that gives Bonferroni intervals (at level \eqn{q} for `"BH"`,
#' at \eqn{q / \sum 1/i}{q / sum(1/i)} for `"BY"`), none of which excludes the
#' null value.
#'
#' The intervals themselves are computed by 'emmeans'
#' (`summary(..., adjust = "none", level = ...)`), so degrees of freedom,
#' back-transformation (for example ratios or odds ratios with
#' `type = "response"`), non-zero `null` values and one-sided intervals
#' (`side`) are handled exactly as in 'emmeans'.
#'
#' Families are defined as in 'emmeans': each group of the `by` variables is a
#' separate family, both for the p-value adjustment and for the interval
#' adjustment. Use `by = NULL` to treat all contrasts as a single family.
#'
#' @section What is and is not guaranteed:
#' The agreement between intervals and adjusted p-values holds by
#' construction (up to floating-point ties at the threshold). Benjamini and
#' Yekutieli (2005) prove that the `"BH"` version controls the FCR, the
#' expected proportion of non-covering intervals among the *selected*
#' (significant) parameters, at level \eqn{q} when the test statistics are
#' independent, and that the `"BY"` version controls it under any dependence
#' structure. Contrasts estimated from the same model are usually correlated,
#' so `method = "BY"` is the conservative choice when a formal guarantee is
#' required. Intervals are also reported for the non-selected contrasts so that
#' the table is complete, but the FCR criterion makes no coverage statement
#' about them.
#'
#' Equivalence tests (`delta`) and cross-family adjustments (`cross.adjust`)
#' are not supported, whether they are passed as arguments or stored in
#' `object` with `update()`.
#'
#' @param object An `emmGrid` object, usually the result of
#'   [emmeans::contrast()] or `pairs()`. For `emm_list` objects (for example
#'   from `emmeans(fit, pairwise ~ x)`), the `contrasts` component is used by
#'   default.
#' @param level Confidence level. The FDR and FCR level is `q = 1 - level`. If
#'   `NULL` (the default), the level stored in `object` is used, which is 0.95
#'   unless it was changed with `update()`.
#' @param method `"BH"` (Benjamini-Hochberg, identical to `adjust = "fdr"` in
#'   'emmeans') or `"BY"` (Benjamini-Yekutieli). If `NULL` (the default), the
#'   method is taken from the `adjust` setting stored in `object`, which must
#'   then be `"fdr"`, `"BH"` or `"BY"`.
#' @param which For `emm_list` objects, the name or index of the component to
#'   use. Defaults to `"contrasts"` when that component exists.
#' @param ... Further arguments passed to [emmeans::summary.emmGrid()], such as
#'   `type = "response"`, `by`, `side`, `null` or `df`. The arguments
#'   `adjust`, `cross.adjust`, `infer` and `delta` (or abbreviations of them)
#'   are not allowed because they would change the family or break the
#'   agreement between intervals and p-values.
#'
#' @return A `summary_emm` data frame, printed in the usual 'emmeans' style,
#'   with FCR-adjusted confidence limits and FDR-adjusted p-values. The number
#'   of tests, the number of rejections and the confidence level used in each
#'   family are stored as an attribute and can be retrieved with [fcr_info()].
#'
#' @references
#' Benjamini, Y. and Yekutieli, D. (2005). False discovery rate-adjusted
#' multiple confidence intervals for selected parameters. *Journal of the
#' American Statistical Association*, 100(469), 71-81.
#' \doi{10.1198/016214504000001907}
#'
#' Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. *Journal of
#' the Royal Statistical Society: Series B*, 57(1), 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' @seealso [fcr_level()], [fcr_info()], [emmeans::summary.emmGrid()].
#'
#' @examplesIf requireNamespace("emmeans", quietly = TRUE)
#' library(emmeans)
#'
#' fit <- lm(weight ~ feed, data = chickwts)
#' prs <- pairs(emmeans(fit, ~ feed), adjust = "fdr")
#'
#' # FDR-adjusted p-values, but Bonferroni-adjusted intervals
#' confint(prs)
#'
#' # FDR-adjusted p-values and FCR-adjusted intervals that agree
#' confint_fcr(prs)
#'
#' # Ratios on the response scale, one family per level of 'percent'
#' fit2 <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' prs2 <- pairs(emmeans(fit2, ~ source | percent), adjust = "fdr")
#' res <- confint_fcr(prs2, type = "response")
#' res
#' fcr_info(res)
#' @export
confint_fcr.emmGrid <- function(object, level = NULL, method = NULL, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for 'emmGrid' objects.",
         call. = FALSE)
  }
  check_dots(names(list(...)))
  check_stored_settings(object@misc)

  if (is.null(level)) {
    level <- object@misc$level
    if (!is.numeric(level) || length(level) != 1L) level <- 0.95
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  q <- 1 - level
  method <- resolve_fcr_method(method, object@misc$adjust)

  # Unadjusted tests define the families and the FDR selection
  tests <- summary(object, infer = c(FALSE, TRUE), adjust = "none",
                   cross.adjust = "none", delta = 0, ...)
  by_vars <- attr(tests, "by.vars")
  family <- family_labels(tests, by_vars)
  rows_by_family <- split(seq_len(nrow(tests)), family)

  p_adj <- rep(NA_real_, nrow(tests))
  info <- lapply(names(rows_by_family), function(fam) {
    rows <- rows_by_family[[fam]]
    p <- tests$p.value[rows]
    if (all(is.na(p))) {
      return(data.frame(family = fam, m = 0L, R = 0L, level = NA_real_))
    }
    fit <- fcr_family(p, q, method)
    data.frame(family = fam, m = fit$m, R = fit$R, level = fit$level)
  })
  info <- do.call(rbind, info)
  for (fam in names(rows_by_family)) {
    rows <- rows_by_family[[fam]]
    p_adj[rows] <- stats::p.adjust(tests$p.value[rows], method = method)
  }

  # Intervals from 'emmeans' at each family's FCR level. Warnings were
  # already shown by the call above, so they are not repeated here.
  out <- NULL
  for (lev in unique(stats::na.omit(info$level))) {
    current <- suppressWarnings(
      summary(object, infer = c(TRUE, TRUE), adjust = "none",
              cross.adjust = "none", delta = 0, level = lev, ...)
    )
    if (is.null(out)) {
      out <- current
      cl_names <- attr(out, "clNames")
    }
    rows <- unlist(rows_by_family[info$family[!is.na(info$level) &
                                                info$level == lev]],
                   use.names = FALSE)
    out[rows, cl_names] <- current[rows, cl_names]
  }
  if (is.null(out)) {
    stop("No non-missing p-values: FCR-adjusted intervals cannot be computed.",
         call. = FALSE)
  }
  no_level <- unlist(rows_by_family[info$family[is.na(info$level)]],
                     use.names = FALSE)
  if (length(no_level) > 0L) {
    out[no_level, cl_names] <- NA_real_
  }

  out$p.value <- p_adj
  attr(out, "adjust") <- if (method == "BH") "fdr" else "BY"
  attr(out, "mesg") <- fcr_messages(attr(out, "mesg"), info, method, q,
                                    by_vars)
  rownames(info) <- NULL
  attr(out, "fcr") <- info
  out
}

#' @rdname confint_fcr.emmGrid
#' @export
confint_fcr.emm_list <- function(object, level = NULL, method = NULL, ...,
                                 which = NULL) {
  if (is.null(which)) {
    if ("contrasts" %in% names(object)) {
      which <- "contrasts"
    } else if (length(object) == 1L) {
      which <- 1L
    } else {
      stop("'object' has several components and none is named ",
           "'contrasts'; choose one with 'which'.", call. = FALSE)
    }
  }
  confint_fcr(object[[which]], level = level, method = method, ...)
}

#' Details of the FCR adjustment
#'
#' Extracts, for each family of tests, the number of tests (`m`), the number
#' of rejected hypotheses (`R`) and the confidence level used for the
#' FCR-adjusted intervals.
#'
#' @param x A result of [confint_fcr()].
#'
#' @return A data frame with columns `family`, `m`, `R` and `level`.
#'
#' @examplesIf requireNamespace("emmeans", quietly = TRUE)
#' library(emmeans)
#' fit <- lm(weight ~ feed, data = chickwts)
#' res <- confint_fcr(pairs(emmeans(fit, ~ feed), adjust = "fdr"))
#' fcr_info(res)
#' @export
fcr_info <- function(x) {
  info <- attr(x, "fcr")
  if (is.null(info)) {
    stop("'x' does not contain FCR information; ",
         "it must be a result of confint_fcr().", call. = FALSE)
  }
  info
}

# ---- Internal helpers -------------------------------------------------------

check_dots <- function(dot_names) {
  forbidden <- c("adjust", "cross.adjust", "infer", "delta")
  dot_names <- dot_names[!is.na(dot_names) & nzchar(dot_names)]
  if (length(dot_names) == 0L) return(invisible(NULL))
  hits <- vapply(dot_names, function(nm) {
    any(startsWith(forbidden, nm))
  }, logical(1))
  if (any(hits)) {
    stop("Argument(s) not allowed in confint_fcr(): ",
         paste0("'", dot_names[hits], "'", collapse = ", "), ". ",
         "The FDR method is set with 'method' (or taken from 'object'); ",
         "equivalence tests and cross-family adjustments are not supported.",
         call. = FALSE)
  }
  invisible(NULL)
}

check_stored_settings <- function(misc) {
  delta <- misc$delta
  if (!is.null(delta) && any(delta != 0)) {
    stop("'object' is set up for equivalence tests (delta != 0), ",
         "which confint_fcr() does not support. ",
         "Use update(object, delta = 0) first.", call. = FALSE)
  }
  cross <- misc$cross.adjust
  if (!is.null(cross) && !identical(tolower(cross), "none")) {
    stop("'object' uses cross.adjust = \"", cross, "\", which ",
         "confint_fcr() does not support. Use update(object, ",
         "cross.adjust = \"none\") first, or by = NULL to pool all ",
         "contrasts into one family.", call. = FALSE)
  }
  invisible(NULL)
}

resolve_fcr_method <- function(method, stored_adjust) {
  if (!is.null(method)) {
    if (!is.character(method) || length(method) != 1L) {
      stop("'method' must be \"BH\" or \"BY\".", call. = FALSE)
    }
    method <- toupper(method)
    if (method == "FDR") method <- "BH"
    if (!method %in% c("BH", "BY")) {
      stop("'method' must be \"BH\" or \"BY\".", call. = FALSE)
    }
    return(method)
  }
  stored <- if (is.null(stored_adjust)) "none" else stored_adjust
  matched <- c("fdr", "bh", "by")[pmatch(tolower(stored),
                                         c("fdr", "bh", "by"))]
  switch(if (is.na(matched)) "other" else matched,
         fdr = "BH",
         bh = "BH",
         by = "BY",
         stop("'object' uses adjust = \"", stored, "\". FCR-adjusted ",
              "intervals correspond to FDR adjustments: create the contrasts ",
              "with adjust = \"fdr\" (or \"BY\"), or set method = \"BH\" ",
              "or \"BY\" in confint_fcr().", call. = FALSE))
}

family_labels <- function(tests, by_vars) {
  if (length(by_vars) == 0L) {
    return(factor(rep("all", nrow(tests))))
  }
  df <- as.data.frame(tests)[by_vars]
  labels <- do.call(paste, c(Map(function(v, x) paste(v, "=", x),
                                 by_vars, df), sep = ", "))
  factor(labels, levels = unique(labels))
}

fcr_messages <- function(mesg, info, method, q, by_vars) {
  pattern <- "^(Confidence level used|P value adjustment|Conf-level adjustment)"
  mesg <- mesg[!grepl(pattern, mesg)]

  sizes <- unique(info$m[info$m > 0])
  size_txt <- if (length(sizes) == 1L) {
    sizes
  } else {
    paste0(min(sizes), "-", max(sizes))
  }
  p_name <- if (method == "BH") "fdr (Benjamini-Hochberg)" else
    "BY (Benjamini-Yekutieli)"
  p_line <- paste0("P value adjustment: ", p_name, " method for ", size_txt,
                   " tests")
  ci_line <- paste0("Conf-level adjustment: FCR (Benjamini & Yekutieli, ",
                    "2005), q = ", format(q))
  fmt_level <- function(i) {
    if (is.na(info$level[i])) return("NA")
    txt <- format(signif(info$level[i], 4))
    if (info$R[i] == 0L) txt <- paste0(txt, " (no rejections: R = 1 used)")
    txt
  }
  if (length(by_vars) == 0L) {
    lev_lines <- paste0("Confidence level used: ", fmt_level(1L))
  } else {
    lev_lines <- c("Confidence level used per family:",
                   vapply(seq_len(nrow(info)), function(i) {
                     paste0("  ", info$family[i], ": ", fmt_level(i))
                   }, character(1)))
  }
  c(lev_lines, ci_line, p_line, mesg)
}
