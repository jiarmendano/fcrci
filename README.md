# fcrci

<!-- badges: start -->
<!-- badges: end -->

**FDR-consistent confidence intervals for multiple comparisons.**

When p-values are adjusted for the false discovery rate (FDR), the confidence
intervals reported next to them are usually adjusted with a different method
(in 'emmeans', `confint()` switches FDR adjustments to Bonferroni). A contrast
can then be significant according to its p-value while its interval includes
the null value.

`fcrci` computes the FCR-adjusted intervals of Benjamini and Yekutieli (2005),
the interval counterpart of the Benjamini-Hochberg procedure. Within each family
of tests, an interval excludes the null value **if and only if** the FDR-adjusted
p-value is at most `1 - level`.

## Installation

The package is not yet on CRAN. To install it from a local copy of the
source:

``` r
# install.packages("remotes")
remotes::install_local("path/to/fcrci", build_vignettes = TRUE)
```

The package itself depends only on base R; 'emmeans' is needed to use
`confint_fcr()` on `emmGrid` objects.

## Example

``` r
library(emmeans)
library(fcrci)

fit <- lm(weight ~ feed, data = chickwts)
prs <- pairs(emmeans(fit, ~ feed), adjust = "fdr")

confint(prs)       # FDR p-values, but Bonferroni intervals
confint_fcr(prs)   # FDR p-values and FCR intervals that agree
```

```
 contrast              estimate   SE df lower.CL upper.CL t.ratio p.value
 casein - horsebean      163.38 23.5 65   113.28   213.49   6.957  <0.0001
 casein - linseed        104.83 22.4 65    57.06   152.61   4.682  <0.0001
 casein - meatmeal        46.67 22.9 65    -2.17    95.52   2.039  0.0570
 ...
 horsebean - linseed     -58.55 23.5 65  -108.66    -8.44  -2.493  0.0228
 ...

Confidence level used: 0.9633
Conf-level adjustment: FCR (Benjamini & Yekutieli, 2005), q = 0.05
P value adjustment: fdr (Benjamini-Hochberg) method for 15 tests
```

With `confint(prs)`, the interval for `horsebean - linseed` is
(-130.1, 13.0) and includes zero, although its FDR-adjusted p-value is 0.023.

Main functions:

* `confint_fcr()`: FCR-adjusted intervals and FDR-adjusted p-values for
  'emmeans' contrasts (`emmGrid` or `emm_list` objects), with one family per
  `by` group, back-transformed ratios and odds ratios (`type = "response"`),
  and the Benjamini-Yekutieli version for dependent tests (`method = "BY"`).
* `fcr_level()`: the FCR-adjusted confidence level for any vector of p-values.
* `fcr_info()`: number of tests, rejections and level used in each family.

See `vignette("fcrci")` for details.

## Development

fcrci was written by Joaquín Armendano, with the assistance of Claude
Code (Anthropic), under manual review.

## References

Benjamini, Y. and Yekutieli, D. (2005). False discovery rate-adjusted multiple
confidence intervals for selected parameters. *Journal of the American
Statistical Association*, 100(469), 71-81.
<https://doi.org/10.1198/016214504000001907>

Jung, K., Friede, T. and Beissbarth, T. (2011). Reporting FDR analogous
confidence intervals for the log fold change of differentially expressed
genes. *BMC Bioinformatics*, 12, 288.
<https://doi.org/10.1186/1471-2105-12-288>
