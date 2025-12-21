# Changelog

All notable changes to `pretest` will be documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)  
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## [0.1.0] - 2025-12-21

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

### Reference
Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test.
*arXiv:2510.26470*. https://arxiv.org/abs/2510.26470
