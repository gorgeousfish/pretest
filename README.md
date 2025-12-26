# pretest

**Conditional Extrapolation Pre-Test for Difference-in-Differences**

[![Stata 17+](https://img.shields.io/badge/Stata-17%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)](CHANGELOG.md)

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
- Survey data with irregular interview schedules (e.g., `nlswork` with gaps)
- Sample restrictions that eliminate entire time periods for one group

**Solution:** Ensure at least one observation per time-treatment cell, or restrict
analysis to time periods with complete coverage.

### Two-Period Designs

This command **cannot be used** for canonical 2×2 DID designs with only two time
periods. For such cases, use standard DID estimators:

```stata
regress y i.treat##i.post
xtdidregress (y) (d), group(id) time(t)
```

## Installation

### From GitHub

```stata
net install pretest, from("https://raw.githubusercontent.com/gorgeousfish/pretest/main/") replace
```

### Local Installation

```stata
net install pretest, from("/path/to/pretest-stata/") replace
```

## Quick Start

```stata
* Load the included example data (California Proposition 99)
findfile prop99_smoking.dta
use "`r(fn)'", clear

* Set panel structure
xtset state year

* Run pre-test with threshold M = 5
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)
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

## Key Formulas

### Pre-test (Theorem 1)

The pre-test indicator is defined as:

> φ = 𝟙{*Ŝ* <sub>pre </sub> > *M*}

where φ = 0 indicates **PASS** (extrapolation justified) and φ = 1 indicates
**FAIL** (extrapolation rejected). See Section 4.2 for the test definition.

### Average DID Estimate

**Important:** The δ̄̂ reported by this package is **not** the traditional ATT.

The DID estimand at time *t* is defined relative to the treatment time *t* <sub>0 </sub>:

> δ̂<sub>t </sub> = (Ȳ <sub>t,D=1 </sub> − Ȳ <sub>t₀,D=1 </sub>) − (Ȳ <sub>t,D=0 </sub> − Ȳ <sub>t₀,D=0 </sub>)

where Ȳ <sub>t,D=d </sub> denotes the sample mean of outcomes for group *D* = *d* at time *t*.

The average DID estimand across post-treatment periods is:

> δ̄̂ = (1/*T* <sub>post </sub>) × Σ <sub>t=t₀</sub><sup>T </sup> δ̂<sub>t </sub>

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

**Iterative mode (Section 5.1):**

> *I* = δ̄̂ ± {κ · *Ŝ* <sub>pre </sub> + *f*(α, Σ̂) / √*n*}

**Overall mode (Appendix C):**

> *I* <sup>Δ </sup> = δ̄̂ ± {*Ŝ* <sup>Δ </sup><sub>pre </sub> + *f* <sup>Δ </sup>(α, Σ̂<sup>Δ </sup>) / √*n*}

The overall mode has no κ multiplier because cumulative violations directly bound the bias (see Appendix C, Proposition).

### κ Constant (Section 3.2, Proposition 1)

> κ = ((1/*T* <sub>post </sub>) · Σ <sub>t=1 </sub><sup>T <sub>post </sub></sup> *t* <sup>q </sup>)<sup>1/q </sup>

where *q* is the Hölder conjugate of *p*, i.e., 1/*p* + 1/*q* = 1.

## Stored Results

### Scalars

| Result              | Description                                                     |
| :------------------ | :-------------------------------------------------------------- |
| `e(S_pre)`        | Estimated pre-treatment severity                                |
| `e(kappa)`        | Bias bound constant κ (iterative mode)                         |
| `e(phi)`          | Pre-test result (0 = pass, 1 = fail, . = data issue or invalid) |
| `e(data_valid)`   | Data validity indicator                                         |
| `e(pretest_pass)` | Pre-test pass indicator                                         |
| `e(delta_bar)`    | Average DID estimate                                            |
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

## Auxiliary Commands

### `pretest_coverage`

Monte Carlo simulation for validating conditional coverage properties.

```stata
pretest_coverage, threshold(0.5) nsims(1000) nobs(200)
```

## References

Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test: Valid Inference When Testing Violations of Parallel Trends for Difference-in-Differences. *arXiv preprint arXiv:2510.26470*. Available at: https://arxiv.org/abs/2510.26470

Rambachan, A., & Roth, J. (2023). A More Credible Approach to Parallel Trends. *Review of Economic Studies*, 90(5), 2555–2591. https://doi.org/10.1093/restud/rdad018

Roth, J. (2022). Pretest with Caution: Event-Study Estimates after Testing for Parallel Trends. *American Economic Review: Insights*, 4(3), 305–322. https://doi.org/10.1257/aeri.20210236

## Authors

**Stata Implementation:**

- Xuanyu Cai
- Wenli Xu

**Methodology:**

- Jonas M. Mikhaeil (Columbia University)
- Christopher Harshaw (Columbia University)

## License

AGPL-3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite both the methodology paper
and the Stata implementation:

```bibtex
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
