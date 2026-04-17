# Changelog

All notable changes to `pretest` will be documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## [0.1.2] - 2026-04-17

### Changed

- **`_pretest_output.ado`**: display labels converted to pure ASCII (e.g. `Alpha:`, `Kappa:`, `Severity (S_pre)`, `delta_bar (vs t0)`, `phi = 1{S_pre > M}`, `CI = delta_bar +/- {kappa*S_pre + f(a,S)/sqrt(n)}`) to avoid SMCL rendering issues on non-UTF8 locales and in log files; `t0` now shows the calendar value rather than the internal 1-based index
- **`_pretest_output.ado`**: `S_pre` precision tightened from 4 decimals to 3 decimals in the Pre-Test Results block to match the scale of the statistic reported in the paper stlog

### Added

- **`_pretest_output.ado`**: new `SE(S_pre)` display line below `Severity (S_pre)` reporting the delta-method standard error of the severity estimate (value passed through from `pretest.ado` via the new `sdpre()` option)
- **`pretest.ado`**: pipes the previously stored `e(S_pre_se)` through to `_pretest_output` via a new internal `sdpre()` option (no public syntax change)

### Notes

- No changes to estimators, interval constructions, or stored results; all `e()` scalars and macros continue to match v0.1.1 numerically
- Accompanying Stata Journal paper aligned with package output: Theorem 2 coverage target stated as the onset-relative DID summary delta_bar, with a remark connecting to the ATT under the onset normalization

---

## [0.1.1] - 2026-04-02

### Fixed

- **`_pretest_graph.ado`**: `saving("file.png", replace)` syntax now correctly parses the filename and replace flag; previously produced a malformed `.gph` file with `, replace` embedded in the filename (BUG-017)
- **`_pretest_output.ado`**: footer note "PASS means φ=0 (extrapolation justified)" is now suppressed when `phi` is missing (data-invalid / no-decision path); replaced with "pre-test decision unavailable (data issue)" (BUG-019)
- **`examples/example_prop99.log`**: regenerated to include all four analyses, remove stale "Coverage check passed" text, and eliminate graph-generation error messages (BUG-016)

---

## [0.1.0] - 2025-12-22

### Added

**Core functionality**

- DID estimators for δ̂ₜ and ν̂ₜ (Section 2.1)
- Severity measurement with *p*-norm, *p* ∈ [1, ∞] (Section 3.1)
- κ constant with closed-form solutions (Section 3.2, Proposition 1)
- ψ function for Monte Carlo critical values (Section 5.1)
- Conditionally valid confidence intervals (Section 5.1, Theorem 2)

**Stata interface**

- Complete command syntax with required and optional arguments
- Formatted output display following Stata conventions
- Comprehensive return values in `e()`
- Event study visualization
- SMCL help documentation (`help pretest`)

**Two violation modes**

- Iterative (default): CI = δ̄̂ ± {κ · Ŝₚᵣₑ + f(α, Σ̂) / √n}
- Overall (Appendix C): CI = δ̄̂ ± {Ŝᐩₚᵣₑ + fᐩ(α, Σ̂ᐩ) / √n}

**Additional features**

- Cluster-robust standard errors
- Coverage simulation command (`pretest_coverage`)

### Requirements

- Stata 17.0 or later
