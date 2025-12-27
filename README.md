# pretest

**Conditional Extrapolation Pre-Test for Difference-in-Differences**

[![Stata 17+](https://img.shields.io/badge/Stata-17%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)](CHANGELOG.md)

![1766840162916](image/README/1766840162916.png)

## Overview

`pretest` implements the conditional extrapolation pre-test framework for
difference-in-differences (DID) designs proposed by Mikhaeil and Harshaw (2025).

Conventional pre-tests for parallel trends suffer from a fundamental limitation:
failing to reject the null hypothesis is uninformative, as the test may simply
lack power. Moreover, conditioning on passing such tests can lead to severe
inference distortions (Roth, 2022). This package provides a principled solution.

Under the **conditional extrapolation assumption**, if pre-treatment violations
fall below an acceptable threshold *M*, extrapolation to post-treatment is justified:

> **Assumption 3 (Conditional Extrapolation):** If *S* <sub>pre </sub> ≤ *M*, then *S* <sub>post </sub> ≤ *S* <sub>pre </sub>.

The package provides:

- An **asymptotically consistent pre-test** for the extrapolation condition
- **Conditionally valid confidence intervals** for the Average Treatment Effect on the Treated (ATT) with guaranteed asymptotic coverage

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
* Load the included example data (California Proposition 99)
webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"
webuse prop99_smoking.dta, clear

* Set panel structure
xtset state year

* Run pre-test with threshold M = 5
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)

* Or use overall mode (less conservative)
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5) overall
```

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
| `treat_time(#)`    | auto    | Treatment time*t* <sub>0 </sub> |
| `p(#)`             | 2       | Severity norm*p* ≥ 1                   |
| `alpha(#)`         | 0.05    | Significance level                        |
| `level(#)`         | 95      | Confidence level (%)                      |
| `cluster(varname)` | —      | Cluster variable for robust SE            |
| `overall`          | off     | Use overall violations mode               |
| `nograph`          | off     | Suppress event study graph                |
| `simulate(#)`      | 5000    | Monte Carlo simulations                   |
| `seed(#)`          | 12345   | Random seed                               |
| `diagnose`         | off     | Display detailed diagnostic information   |

### Graph Customization

| Option                          | Description                           |
| :------------------------------ | :------------------------------------ |
| `ci_opt_pass(string)`         | Override CI style when pretest passes |
| `ci_opt_fail(string)`         | Override CI style when pretest fails  |
| `line_opt_m(string)`          | Override threshold M line style       |
| `marker_opt_pre(string)`      | Override pre-treatment marker style   |
| `marker_opt_post(string)`     | Override post-treatment marker style  |
| `scheme()`, `title()`, etc. | Any standard Stata twoway_options     |

**Note:** Element-specific options (e.g., `ci_opt_pass`) completely replace default styling when specified.

## Key Formulas

### Pre-test (Theorem 1)

The pre-test indicator is defined as:

> φ = 𝟙{*Ŝ* <sub>pre </sub> > *M*}

where φ = 0 indicates **PASS** (extrapolation justified) and φ = 1 indicates
**FAIL** (extrapolation rejected).

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

1. The CI bounds account for potential bias via κ · Ŝ_pre
2. The interpretation is: "treatment effect relative to treatment onset"

For traditional ATT comparison, use `e(ci_conv_lower)` and `e(ci_conv_upper)`.

### Conditional Confidence Interval (Theorem 2)

 **1. Iterative mode (Default):**

> *I* = δ̄̂ ± {κ · *Ŝ* <sub>pre </sub> + *f*(α, Σ̂) / √*n*}

 Bias bound includes the multiplier κ ≥ 1.

 **2. Overall mode:**

> *I*<sup>Δ</sup> = δ̄̂ ± {*Ŝ*<sup>Δ</sup><sub>pre</sub> + *f*<sup>Δ</sup>(α, Σ̂<sup>Δ</sup>) / √*n*}

 Bias bound uses *no multiplier* (κ = 1).

### κ Constant (Iterative Mode Only)

> κ = ((1/*T*<sub>post</sub>) · Σ<sub>t=1</sub><sup>T<sub>post</sub></sup> *t*<sup>q</sup>)<sup>1/q</sup>

 where *q* is the Hölder conjugate of *p*. κ captures the worst-case accumulation of iterative violations over time.

- For *T* <sub>post </sub> > 1, κ > 1.
- For *p* = 2 and large *T* <sub>post </sub>, κ grows with √*T* <sub>post </sub>.
- **Overall Mode:** κ is not used (effectively κ = 1), yielding narrower intervals.

## Stored Results

### Scalars

| Result              | Description                                                     |
| :------------------ | :-------------------------------------------------------------- |
| `e(S_pre)`        | Estimated pre-treatment severity                                |
| `e(S_pre_se)`     | Standard error of S_pre (Delta method)                          |
| `e(kappa)`        | Bias bound constant κ (iterative mode)                         |
| `e(phi)`          | Pre-test result (0 = pass, 1 = fail, . = data issue or invalid) |
| `e(data_valid)`   | Data validity indicator                                         |
| `e(pretest_pass)` | Pre-test pass indicator                                         |
| `e(delta_bar)`    | Average DID estimate                                            |
| `e(se_delta_bar)` | Standard error of average DID estimate                          |
| `e(ci_lower)`     | Conditional CI lower bound                                      |
| `e(ci_upper)`     | Conditional CI upper bound                                      |
| `e(T)`            | Total time periods                                              |
| `e(T_pre)`        | Pre-treatment periods                                           |
| `e(T_post)`       | Post-treatment periods                                          |
| `e(N)`            | Number of observations                                          |

### Matrices

| Result       | Description                                                  |
| :----------- | :----------------------------------------------------------- |
| `e(nu)`    | Iterative violations (*T* <sub>pre </sub>−1 × 1) |
| `e(delta)` | DID estimates (*T* <sub>post </sub> × 1)          |
| `e(theta)` | Full parameter vector θ̂                                   |
| `e(Sigma)` | Asymptotic covariance matrix                                 |
| `e(b)`     | Coefficient vector                                           |
| `e(V)`     | Variance matrix                                              |

## Mode Selection: Iterative vs. Overall

The package offers two assumptions about parallel trend violations, which have different sensitivities:

| Feature               | Iterative Mode (Default)                                        | Overall Mode (`overall`)                               |
| :-------------------- | :-------------------------------------------------------------- | :------------------------------------------------------- |
| **Assumption**  | Violations accumulate period-to-period                          | Violations are bounded by cumulative total               |
| **Sensitivity** | Sensitive to**volatility/noise** (sharp changes)          | Sensitive to**drift/trend** (long-term divergence) |
| **Blind Spot**  | May pass smooth linear trends (constant small changes)          | May fail even if period-to-period changes are small      |
| **Bias Bound**  | Scaled by κ (proportional to √*T* <sub>post </sub>) | **No multiplier** (κ = 1)                         |
| **CI Width**    | Generally Wider (accounts for worst-case accumulation)          | Generally Narrower (assumes bounded total error)         |

**Recommendation:**

1. **Start with Iterative Mode** as it is the standard, robust approach.
2. **Check Overall Mode if:**
   - You suspect a **linear trend** or long-term drift (Iterative might incorrectly pass).
   - The Iterative results are too conservative (wide CIs) despite visually decent parallel trends.
3. **If Iterative PASSES but Overall FAILS:** This strongly suggests the presence of a smooth linear trend difference between groups. Extrapolation is risky unless you conceptually allow for this drift to continue.

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
xtset id time

* Run pre-test
pretest y, treatment(treat) time(time) threshold(0.5) treat_time(6)
```

## Future Roadmap

The development team is evaluating the addition of **Threshold Sensitivity Analysis** in future versions.
Current methodology treats $M$ as a fixed gatekeeper. Future updates may introduce a sensitivity analysis mode that:

- Visualizes how the Conditional Confidence Interval (CI) varies across a continuous range of $M$ values.
- Adopts a "Partial Identification" perspective to show the robustness of conclusions to different assumptions about the maximum acceptable violation.

## References

Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test: Valid Inference When Testing Violations of Parallel Trends for Difference-in-Differences. *arXiv preprint arXiv:2510.26470*. Available at: https://arxiv.org/abs/2510.26470

Rambachan, A., & Roth, J. (2023). A More Credible Approach to Parallel Trends. *Review of Economic Studies*, 90(5), 2555–2591. https://doi.org/10.1093/restud/rdad018

Roth, J. (2022). Pretest with Caution: Event-Study Estimates after Testing for Parallel Trends. *American Economic Review: Insights*, 4(3), 305–322. https://doi.org/10.1257/aeri.20210236

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of MacauEmail: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
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

> Cai, X., & Xu, W. (2025). *pretest: Stata module to implement the conditional extrapolation pre-test for difference-in-differences* (Version 0.1.0) [Computer software]. GitHub. https://github.com/gorgeousfish/pretest
>
> Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test: Valid Inference when Testing Violations of Parallel Trends for Difference-in-Differences. *arXiv preprint arXiv:2510.26470*. https://arxiv.org/abs/2510.26470

**BibTeX:**

```bibtex
@software{pretest2025stata,
      title={pretest: Stata module to implement the conditional extrapolation pre-test for difference-in-differences},
      author={Xuanyu Cai and Wenli Xu},
      year={2025},
      version={0.1.0},
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
