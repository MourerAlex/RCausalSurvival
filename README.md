# CausalSurvival

> **Status**: pre-implementation. Skeleton only. See `dev/CAUSAL_SURVIVAL_SPEC.md` in the [CausalCompetingRisks](https://github.com/MourerAlex/CausalCompetingRisks) repository for the full specification.

R package for causal inference on single-event survival outcomes, using discrete-time pooled logistic regression. Provides parametric g-formula and inverse probability weighting estimators for cumulative incidence under static, baseline-only treatment regimes, with bootstrap confidence intervals and identifying-assumption accessors.

Designed as the foundation of a two-package ecosystem: [CausalCompetingRisks](https://github.com/MourerAlex/CausalCompetingRisks) extends to competing events via the separable-effects framework.

## Installation

Not yet on CRAN. Development version:

```r
# install.packages("remotes")
remotes::install_github("MourerAlex/RCausalSurvival")
```

## Usage

```r
library(CausalSurvival)

# --- Synthetic data ---------------------------------------------------------
set.seed(42); n <- 200
df <- data.frame(id = seq_len(n), L1 = rnorm(n), L2 = rbinom(n, 1, 0.4))
df$A      <- rbinom(n, 1, plogis(-0.2 + 0.3 * df$L1 + 0.5 * df$L2))
df$lambda <- plogis(-2 + 0.3 * df$A + 0.2 * df$L1 + 0.4 * df$L2)
df$time   <- pmin(rgeom(n, df$lambda) + 1L, 10L)
df$status <- as.integer(df$time < 10L)

pt <- to_person_time(df, id = "id", time = "time", status = "status",
                     treatment = "A", covariates = c("L1", "L2"),
                     cut_points = 20)

# --- g-formula --------------------------------------------------------------
fit_g <- causal_survival(pt, method = "gformula")
print(fit_g)
summary(fit_g)
print(causal_assumptions(fit_g))

# --- IPW, default km estimator ----------------------------------------------
fit_i <- causal_survival(pt, method = "ipw", truncate = c(0.01, 0.99))
print(fit_i)

# --- IPW, msm estimator (arg renamed: .ipw_engine -> .ipw_estimator) --------
fit_m <- causal_survival(pt, method = "ipw", truncate = c(0.01, 0.99),
                         .ipw_estimator = "msm")
print(fit_m)

# --- Accessors --------------------------------------------------------------
print(causal_risk(fit_g, "incidence"))
print(causal_contrast(fit_g))  # emits the loud ci = NULL warning

# --- Bootstrap + contrast with CI -------------------------------------------
boot_g <- bootstrap(fit_g, n_boot = 500, alpha = 0.05, seed = 1)
print(boot_g)
print(causal_contrast(fit_g, ci = boot_g))
summary(fit_g, ci = boot_g)

# --- Risk-table accessor ----------------------------------------------------
print(causal_risk_table(fit_g, count = "at_risk"))

# --- Plot -------------------------------------------------------------------
# Each fit needs its OWN bootstrap. The replicates are method-specific —
# bootstrap(fit_g) stores g-formula CIFs, so using boot_g with fit_i's plot
# would silently show g-formula bands on an IPW curve.
# plot.1 bootstrap for each fit
boot_i <- bootstrap(fit_i, n_boot = 500, alpha = 0.05, seed = 1)
boot_m <- bootstrap(fit_m, n_boot = 500, alpha = 0.05, seed = 1)

# plot.2 different fits plots with different options
plot(causal_risk(fit_g, "incidence", ci = boot_g), risk_table = "at_risk")
plot(causal_risk(fit_i, "survival", ci = boot_i), risk_table = "events_y")
plot(causal_risk(fit_m, "survival", ci = boot_m))

# plot.3 contrast plot
plot(causal_contrast(fit_g, ci = boot_g))

# plot.4 fit plot with stacked tables and tables only
plot(causal_risk(fit_g),
     risk_table = c("at_risk", "events_y", "censored"))
plot(causal_risk(fit_g),
     risk_table = c("at_risk", "events_y", "censored"),
     curves     = FALSE)
```

See `vignette("getting-started")` once available.

## Notes & gotchas

Short answers to the questions that come up most. A full vignette is planned.

**Two kinds of censoring — which is which?**
`to_person_time()` splits censoring into two mechanisms, labelled per subject
via the `ipcw` argument (a logical vector; `TRUE` = the first kind):

- `cond_indep_cens` — *conditionally-independent* censoring: independent of
  the outcome **only after** conditioning on covariates `L` and treatment `A`
  (e.g. loss-to-follow-up, treatment switching). It is modelled and corrected
  by IPCW.
- `indep_cens` — *independent* censoring: assumed independent of the outcome
  **without** conditioning (administrative end-of-study is the usual example).
  It gets weight 1 (no model).

Label by the *independence assumption you are willing to make*, not by the
administrative reason — administrative censoring is the typical example of
`indep_cens` but is not automatically independent.

**Why is `w_cens` ≈ 1 even though a censoring model was fit?**
With `ipcw = TRUE` the package fits the `cond_indep_cens` hazard model. If your
data have no conditionally-independent censoring inside `(0, T_max]` (e.g. every
censored subject is administratively truncated at `T_max`), that model sees no
events, predicts ~0 hazard, and the IPCW collapses to 1. Not a bug — it reflects
that there was nothing for IPCW to correct.

**What does `cut_points` do?**
Discretizes follow-up `(0, T_max]` into intervals. `NULL` → 12 equal-width
intervals; a single integer `N` → `N` equal-width intervals; a numeric vector →
explicit interior cut points. `k` is the integer interval index (1, 2, …), not
elapsed time.

**Time-varying treatment / covariates?**
Not in v1. Treatment is **point (baseline)**: `A_0` is carried across all
intervals, and the propensity model is fit on the baseline row only. The KM IPW
engine cannot support time-varying `A` (the arm risk set is undefined); the MSM
engine is the seam for a future extension.

**Is it pooled logistic, or a cubic in time?**
Both. The hazard models are pooled logistic regressions (one logistic fit over
all person-time rows); time enters flexibly as `k + I(k^2) + I(k^3)` rather than
as interval dummies. Default: `<event> ~ A + k + I(k^2) + I(k^3) + covariates`.
Override any model with `formulas = list(y = ..., c = ..., A = ...)`.

**`fit$models$A` vs `fit$models$A_num`?**
The IPW treatment weight is **stabilized**: `A` is the denominator propensity
model `A ~ L`, `A_num` is the numerator marginal model `A ~ 1`. The weight is
`w_a = P(A_obs | A~1) / P(A_obs | A~L)`, then truncated when `truncate` is set.
Inspect raw vs truncated weights in `fit$weights$weight_summary`.

**How does the g-formula handle censoring?**
It fits no censoring model. The Y-hazard is fit on uncensored at-risk rows and
the recursion targets the `c = 0` counterfactual under **independent censoring
given `L_0`**. The `cond_indep_cens` / `indep_cens` split only matters for IPW.

**Is bootstrap the only way to get confidence intervals?**
In v1, yes (analytic / influence-function CIs are planned). `bootstrap()`
returns a standalone object rather than baking CIs into the fit so that: the
fit stays pure point estimates; bootstrapping (optional and slow) is opt-in with
its own `n_boot` / `alpha` / `seed`; and one bootstrap is reused across
`causal_risk()`, `causal_contrast()`, `plot()`, and `summary()`. Replicates are
method-specific, so **each fit needs its own bootstrap** (a g-formula bootstrap
on an IPW curve would silently mislabel the bands).

## License

MIT © Alex Mourer
