# pretest

**Thresholded DID Diagnostics for Conditional Extrapolation**

[![Stata 17+](https://img.shields.io/badge/Stata-17%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.1](https://img.shields.io/badge/Version-0.1.1-green.svg)](CHANGELOG.md)

![1766840162916](image/README/1766840162916.png)

## Overview

`pretest` asks whether the observed pre-treatment path is mild enough to
carry `delta_bar` into conditional extrapolation in common-treatment-time
block-adoption designs. One run reports pre-treatment severity `S_pre`, a
pass indicator, `delta_bar`, and, only when `pretest_pass = 1`, pass-case
interval endpoints read as centered on `delta_bar` and informative about the
average ATT only under the conditional extrapolation assumption.

The shipped Prop99 example shows the decision margin. Lowering the threshold
from 5 to 1 leaves `delta_bar = -14.4450` and `S_pre = 2.1959` unchanged but
moves `pretest_pass` from 1 to 0. The threshold changes what the user may
claim before it changes the DID summary.

In practice, the command separates three choices that applied DID work often
blurs: threshold choice, sample choice, and violation-summary choice.

Read a run in four steps: severity first, gate second, `delta_bar` third, and
pass-case interval endpoints last. When those endpoints appear, read them as
centered on `delta_bar`. Treat them as informative about the average ATT only
under the conditional extrapolation assumption and the asymptotic conditions
in the method paper.

Users choose the threshold under the **conditional extrapolation assumption**:

> **Assumption 3 (Conditional Extrapolation):** If *S* <sub>pre </sub> ≤ *M*, then *S* <sub>post </sub> ≤ *S* <sub>pre </sub>.

A typical run gives users:

- The **pre-treatment severity** `S_pre`
- A **pass indicator** for the chosen threshold
- The DID summary **`delta_bar`**
- **Pass-case interval endpoints** for that conditional reading rule

### First Run

```stata
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)
```

Use `overall` only when the cumulative-violation path is the empirical
question.

## Requirements

**Before using this package, ensure your data meets the following requirements:**

| Requirement                         | Description                                                                                                                                                                 |
| :---------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Minimum 3 time periods**    | *T* <sub>pre </sub> ≥ 2. At least two pre-treatment periods are required because iterative violations ν̂<sub>t </sub> are only defined for *t* ≥ 2. |
| **Block adoption design**     | All treated units must receive treatment at the same time*t* <sub>0 </sub>. Staggered adoption designs are **not** supported.                               |
| **Binary treatment**          | Treatment indicator must be coded as 0 (control) or 1 (treated).                                                                                                            |
| **Complete time-group cells** | Each time period must contain observations in both treatment and control groups.                                                                                            |

### Data Completeness

When some time periods lack observations for either group, the covariance matrix
cannot be computed. In such cases:

- `e(phi)` = . (missing, indicating data issue)
- `e(data_valid)` = 0
- `e(S_pre)`, `e(f_alpha)`, `e(ci_lower)`, `e(ci_upper)` will be missing (.)

**Common causes:**

- Missing values in treatment or outcome variables creating empty cells
- Survey data with irregular interview schedules
- Sample restrictions that eliminate entire time periods for one group

**Solution:** Ensure at least one observation per time-treatment cell, or restrict
analysis to time periods with complete coverage.

### Two-Period Designs

This command **cannot be used** for canonical 2×2 DID designs with only two time
periods.

## Installation

### From GitHub

```stata
net install pretest, from("https://raw.githubusercontent.com/gorgeousfish/pretest/main") replace
```

### From SSC (Coming Soon)

This package will be available on the SSC archive in the near future.

## Quick Start

```stata
* Load the example data (prefer a local shipped or net-get copy)
capture noisily findfile prop99_smoking.dta
if _rc == 0 {
    use "`r(fn)'", clear
}
else {
    webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"
    webuse prop99_smoking.dta, clear
}

* Set panel structure
xtset state year

* Run pre-test with threshold M = 5
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)

* Use overall only when the cumulative-violation path is the question
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5) overall
```

**Warning for clustered runs:** Before clustered runs, subset the analysis
sample explicitly and rerun `pretest` on that fixed sample; inline `[if] [in]`
restrictions are not a supported substitute for this workflow.

`net get pretest` downloads loose working-directory copies of
`prop99_smoking.dta` and `example_prop99.do`. A clean `net install pretest`
alone does not place discoverable example assets on the installed adopath, so
the GitHub `webuse` path is a fallback rather than the canonical local route.

## Syntax

```stata
pretest depvar , treatment(varname) time(varname) threshold(#) [options]
```

### Required

| Option                 | Description                             |
| :--------------------- | :-------------------------------------- |
| `treatment(varname)` | Binary treatment indicator (0/1)        |
| `time(varname)`      | Time variable                           |
| `threshold(#)`       | Acceptable violation threshold*M* > 0 |

### Optional

| Option               | Default | Description                               |
| :------------------- | :------ | :---------------------------------------- |
| `treat_time(#)`    | auto    | Common treatment-onset time*t* <sub>0 </sub>; auto-detected only from a common 0→1 switch path |
| `p(#)`             | 2       | Severity norm*p* ≥ 1                   |
| `alpha(#)`         | `1-c(level)/100` | Nominal significance level for interval construction |
| `level(#)`         | `c(level)` | Nominal confidence level (%) for interval construction |
| `cluster(varname)` | —      | Cluster variable for robust SE after explicit sample subsetting and missing-row cleanup |
| `overall`          | off     | Use cumulative-violation path             |
| `nograph`          | off     | Suppress event study graph                |
| `simulate(#)`      | 5000    | Monte Carlo simulations                   |
| `seed(#)`          | 12345   | Random seed                               |
| `diagnose`         | off     | Display detailed diagnostic information   |

### Graph Customization

| Option                          | Description                           |
| :------------------------------ | :------------------------------------ |
| `ci_opt_pass(string)`         | Override pre-treatment CIs and pass-case post-treatment CIs |
| `ci_opt_fail(string)`         | Override CI style in non-pass cases   |
| `line_opt_m(string)`          | Override threshold M line style       |
| `marker_opt_pre(string)`      | Override pre-treatment marker style   |
| `marker_opt_post(string)`     | Override post-treatment marker style  |
| `scheme()`, `title()`, etc. | Any standard Stata twoway_options     |

**Note:** Element-specific options replace the default styling for the layers
they control. In particular, `ci_opt_pass()` controls the pre-treatment CI
layer and, when the run passes, the post-treatment pass-case endpoints.

If `alpha()` and `level()` are both omitted, `pretest` inherits Stata's current
confidence level. If both are supplied, `alpha()` takes precedence.

## Key Formulas

### Pre-test Indicator

At the interface level, use the gate:

> φ = 𝟙{*Ŝ* <sub>pre </sub> > *M*}

Here φ = 0 means the run passes the chosen threshold and φ = 1 means it
fails that threshold. The main reporting flag is `pretest_pass`: it equals 1
when `phi = 0` and 0 in all non-pass cases. When `phi` is missing because
the data are invalid, `pretest_pass` still remains 0, so read 0 as non-pass
and use `phi` together with `data_valid` to separate a threshold fail from an
invalid-data run.

### Average DID Estimate

**Important:** The δ̄̂ reported by this package is **not** the traditional ATT.

The DID estimand at time *t* is defined relative to the treatment time *t* <sub>0 </sub>:

> δ̂<sub>t </sub> = (Ȳ <sub>t,D=1 </sub> − Ȳ <sub>t₀,D=1 </sub>) − (Ȳ <sub>t,D=0 </sub> − Ȳ <sub>t₀,D=0 </sub>)

where Ȳ <sub>t,D=d </sub> denotes the sample mean of outcomes for group *D* = *d* at time *t*.

The average DID estimand across post-treatment periods is:

> δ̄̂ = (1/*T*<sub>post</sub>) × Σ<sub>t=t₀</sub><sup>T</sup> δ̂<sub>t</sub>

**Key differences from traditional DID:**

| Aspect                     | Paper's δ̄̂                                     | Traditional ATT        |
| :------------------------- | :------------------------------------------------- | :--------------------- |
| Reference point            | Treatment time*t* <sub>0 </sub>          | Pre-treatment average  |
| δ̂<sub>t₀</sub> | Always 0 (by construction)                         | N/A                    |
| Interpretation             | Incremental change from*t* <sub>0 </sub> | Total treatment effect |

**Example:** If treatment effect is constant at 2.0 per period:

- Traditional ATT ≈ 2.0 (total effect)
- Paper's δ̄̂ ≈ 0 (no incremental change after t₀)

**Why this definition?** The paper's δ̄̂ is designed for the conditional extrapolation framework, where:

1. In iterative mode, the bias bound uses κ · Ŝ_pre; in overall mode, the
   command switches to the cumulative-violation object Ŝ^Δ_pre rather than
   simply reusing the iterative path with κ turned off
2. The summary tracks incremental change from treatment onset rather than the
   ATT level

For comparison with a conventional DID interval for `delta_bar` under exact
parallel trends, use `e(ci_conv_lower)` and `e(ci_conv_upper)`; for ATT-level
estimation, use a standard DID estimator.

### Pass-Case Conditional Interval (Theorem 2)

 **1. Iterative mode (Default):**

> *I* = δ̄̂ ± {κ · *Ŝ* <sub>pre </sub> + *f*(α, Σ̂) / √*n*}

 Bias bound includes the multiplier κ ≥ 1.

 **2. Overall mode:**

> *I*<sup>Δ</sup> = δ̄̂ ± {*Ŝ*<sup>Δ</sup><sub>pre</sub> + *f*<sup>Δ</sup>(α, Σ̂<sup>Δ</sup>) / √*n*}

 Bias bound uses the cumulative-violation object *Ŝ*<sup>Δ</sup><sub>pre</sub>, so
 overall mode is a different inferential path rather than iterative mode with
 κ switched off.

### κ Constant (Iterative Mode Only)

> κ = ((1/*T*<sub>post</sub>) · Σ<sub>t=1</sub><sup>T<sub>post</sub></sup> *t*<sup>q</sup>)<sup>1/q</sup>

 where *q* is the Hölder conjugate of *p*. κ captures the worst-case accumulation of iterative violations over time.

- For *T* <sub>post </sub> > 1, κ > 1.
- For *p* = 2 and large *T* <sub>post </sub>, κ grows with √*T* <sub>post </sub>.
- **Overall Mode:** κ is not used (effectively κ = 1), so the bias bound and
  interval construction follow the cumulative-violation path instead.

## Stored Results

### Scalars

| Result              | Description                                                     |
| :------------------ | :-------------------------------------------------------------- |
| `e(S_pre)`        | Estimated pre-treatment severity                                |
| `e(S_pre_se)`     | Standard error of S_pre (Delta method)                          |
| `e(kappa)`        | Bias-bound scalar κ (posted as 1 in overall mode)             |
| `e(phi)`          | Pre-test result (0 = pass, 1 = fail, . = data issue or invalid) |
| `e(data_valid)`   | Data validity indicator                                         |
| `e(pretest_pass)` | Pre-test pass indicator (1 = pass, 0 = non-pass; read with `e(phi)` and `e(data_valid)` when needed) |
| `e(delta_bar)`    | Average DID change relative to `t0`                             |
| `e(se_delta_bar)` | Standard error of average DID estimate                          |
| `e(ci_lower)`     | Pass-case interval lower endpoint; read as part of a `delta_bar`-centered interval when `pretest_pass = 1` and as informative about the average ATT only under the conditional extrapolation assumption |
| `e(ci_upper)`     | Pass-case interval upper endpoint; read as part of a `delta_bar`-centered interval when `pretest_pass = 1` and as informative about the average ATT only under the conditional extrapolation assumption |
| `e(T)`            | Total time periods                                              |
| `e(t0)`           | Treatment time on internal consecutive index                    |
| `e(t0_orig)`      | Treatment time in original calendar units (e.g., the year supplied to `treat_time()`) |
| `e(T_pre)`        | Pre-treatment periods                                           |
| `e(T_post)`       | Post-treatment periods                                          |
| `e(N)`            | Number of observations                                          |

### Matrices

| Result       | Description                                                  |
| :----------- | :----------------------------------------------------------- |
| `e(nu)`    | Pre-treatment violations (*T* <sub>pre </sub>−1 × 1); iterative ν̂_t in iterative mode, cumulative ν̄_t in overall mode |
| `e(delta)` | DID estimates (*T* <sub>post </sub> × 1)          |
| `e(theta)` | Full parameter vector θ̂; pre-treatment block is ν̂_t (iterative mode) or ν̄_t (overall mode) |
| `e(Sigma)` | Asymptotic covariance matrix; in overall mode stored in overall-frame coordinates (A·Σ·A') |

Compatibility returns such as `e(b)` and `e(V)` remain for Stata
postestimation plumbing. They mirror `delta_bar` and its variance, so use
`e(delta_bar)` as the substantive DID summary.

## Mode Selection: Iterative vs. Overall

The package offers two violation summaries for parallel-trend problems, which
have different sensitivities:

| Feature               | Iterative Mode (Default)                                        | Overall Mode (`overall`)                               |
| :-------------------- | :-------------------------------------------------------------- | :------------------------------------------------------- |
| **Assumption**  | Violations accumulate period-to-period                          | Violations are bounded by cumulative total               |
| **Sensitivity** | Sensitive to**volatility/noise** (sharp changes)          | Sensitive to**drift/trend** (long-term divergence) |
| **Blind Spot**  | May pass smooth linear trends (constant small changes)          | May fail even if period-to-period changes are small      |
| **Bias Bound**  | Scaled by κ (proportional to √*T* <sub>post </sub>) | Uses the cumulative-violation path; no iterative κ multiplier |
| **CI Width**    | Depends on the κ-scaled bias bound                             | Depends on the cumulative-violation path                 |

**How to compare the two modes:**

1. Start with iterative mode when period-to-period volatility is the main concern.
2. Check `overall` when cumulative drift is the empirical question.
3. If one mode passes and the other does not, treat that split as evidence that the two violation summaries answer different questions, not as proof that one mode is uniformly more permissive or more informative.

## Example

```stata
* Simulated panel data
clear
set seed 12345
set obs 500
gen id = ceil(_n/10)
gen time = mod(_n-1, 10) + 1
gen treat = (id <= 25)
gen y = rnormal() + treat*(time >= 6)*0.5
```

## Future Roadmap

The development team is evaluating the following extensions for future versions:

- **Triple Difference-in-Differences (DDD)**: Extending the conditional extrapolation framework to DDD designs with an additional grouping dimension.
- **Staggered Adoption Support**: Extending the conditional extrapolation framework to staggered treatment adoption designs via cohort stacking.
- **Covariate Adjustment**: Adding support for control variables in the estimation.
- **Threshold Sensitivity Analysis**: Visualizing how the pass-case conditional interval varies across a continuous range of $M$ values.

## References

Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test: Valid Inference When Testing Violations of Parallel Trends for Difference-in-Differences. *arXiv preprint arXiv:2510.26470*. Available at: https://arxiv.org/abs/2510.26470

Rambachan, A., & Roth, J. (2023). A More Credible Approach to Parallel Trends. *Review of Economic Studies*, 90(5), 2555–2591. https://doi.org/10.1093/restud/rdad018

Roth, J. (2022). Pretest with Caution: Event-Study Estimates after Testing for Parallel Trends. *American Economic Review: Insights*, 4(3), 305–322. https://doi.org/10.1257/aeri.20210236

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **Jonas M. Mikhaeil**, Department of Statistics, Columbia University
- **Christopher Harshaw**, Department of Statistics, Columbia University

## License

AGPL-3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite both the methodology paper
and the Stata implementation:

**APA Format:**

> Cai, X., & Xu, W. (2025). *pretest: Thresholded DID diagnostics for conditional extrapolation* (Version 0.1.1) [Computer software]. GitHub. https://github.com/gorgeousfish/pretest
>
> Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test: Valid Inference when Testing Violations of Parallel Trends for Difference-in-Differences. *arXiv preprint arXiv:2510.26470*. https://arxiv.org/abs/2510.26470

**BibTeX:**

```bibtex
@software{pretest2025stata,
      title={pretest: Thresholded DID diagnostics for conditional extrapolation},
      author={Xuanyu Cai and Wenli Xu},
      year={2025},
      version={0.1.1},
      url={https://github.com/gorgeousfish/pretest}
}

@misc{mikhaeil2025defensepretestvalidinference,
      title={In Defense of the Pre-Test: Valid Inference when Testing Violations 
             of Parallel Trends for Difference-in-Differences}, 
      author={Jonas M. Mikhaeil and Christopher Harshaw},
      year={2025},
      eprint={2510.26470},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2510.26470}
}
```

## Recommended Resources

For beginners in causal inference and econometrics:

- [Causal Inference for the Brave and True](https://matheusfacure.github.io/python-causality-handbook/landing-page.html) - An excellent introductory tutorial on causal inference by Matheus Facure
- [Causal Inference for the Brave and True (Chinese Edition)](https://ci-book.huangwz.com/intro) - Chinese translation by Wenzhe Huang and Wenli Xu
