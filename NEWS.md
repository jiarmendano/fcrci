# fcrci 0.1.0

* Initial release.
* `confint_fcr()` generic with methods for 'emmeans' `emmGrid` and
  `emm_list` objects:
  FCR-adjusted confidence intervals (Benjamini and Yekutieli, 2005) that agree
  with FDR-adjusted p-values, one family per `by` group.
* `fcr_level()` computes the FCR-adjusted confidence level from a vector of
  p-values, for the Benjamini-Hochberg and Benjamini-Yekutieli procedures.
  `method = "fdr"` is accepted as a synonym of `"BH"`.
* `fcr_info()` returns the number of tests, rejections and the level used in
  each family.
