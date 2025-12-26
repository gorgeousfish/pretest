*! _pretest_output.ado
*! Version 0.1.0
*!
*! Title: Formatted Output Display for Conditional Extrapolation Pre-Test
*!
*! Description:
*!   Internal subroutine that displays formatted results for the pretest
*!   command following standard Stata econometric output conventions.
*!
*! Output Structure:
*!   Section 1 - Header: Command title and citation
*!   Section 2 - Model Specification: Threshold M, mode, alpha, p-norm
*!   Section 3 - Time Structure: T, t0, T_pre, T_post, sample size
*!   Section 4 - Pre-Test Results: Severity estimate, kappa, critical value
*!   Section 5 - Treatment Effects: DID estimate and confidence intervals
*!   Section 6 - Footer: Methodological notes and formulas
*!
*! Numeric Formatting Rules:
*!   Standard values:     %9.3f (three decimal places)
*!   Large values (>=1000): %9.0fc (comma-separated integers)
*!   Extreme values:      %9.3e (scientific notation for |x|>=1e6)
*!   Small values:        %9.6f (six decimals for |x|<0.0005)
*!
*! Econometric Background:
*!   The conditional extrapolation assumption (Assumption 3 in the paper)
*!   states: If S_pre <= M, then S_post <= S_pre, where S denotes the
*!   severity of parallel trends violations measured by the p-norm.
*!
*!   Under this assumption, the confidence interval for the ATT is:
*!     Iterative mode: CI = delta_bar +/- {kappa * S_pre + f(alpha,Sigma)/sqrt(n)}
*!     Overall mode:   CI = delta_bar +/- {S_pre + f(alpha,Sigma)/sqrt(n)}
*!
*!   where kappa = (1/T_post * sum_{t=1}^{T_post} t^q)^{1/q} with 1/p + 1/q = 1.
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test:
*!   Valid Inference when Testing Violations of Parallel Trends for
*!   Difference-in-Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

program define _pretest_output, rclass
    version 17.0
    
    // Uses rclass to preserve caller's e() results from estimation commands
    syntax, ///
        THRESHOLD(real) ///
        MODE(string) ///
        ALPHA(real) ///
        LEVEL(real) ///
        T(integer) ///
        t0(integer) ///
        T_pre(integer) ///
        T_post(integer) ///
        N(integer) ///
        P(real) ///
        S_pre(real) ///
        KAPPA(real) ///
        F_alpha(real) ///
        PHI(real) ///
        DELTA_bar(real) ///
        [CI_lower(real -999) ///
         CI_upper(real -999) ///
         CI_conv_lower(real -999) ///
         CI_conv_upper(real -999) ///
         SE_delta(real -999) ///
         SIMS(integer 5000) ///
         SEED(integer 12345)]
    
    // =========================================================================
    // Section 1: Header
    // Display command title and literature reference
    // =========================================================================
    di as text ""
    di as text "{hline 72}"
    di as text "{bf:Conditional Extrapolation Pre-Test for Difference-in-Differences}"
    di as text "Mikhaeil & Harshaw (2025)"
    di as text "{hline 72}"
    
    // =========================================================================
    // Section 2: Model Specification
    // Display user-specified parameters: threshold M, violation mode, alpha
    // =========================================================================
    di as text ""
    di as text "{bf:Model Specification}"
    di as text "{hline 72}"
    
    // Column alignment for parameter table
    local col1 = 4
    local col2 = 28
    local col3 = 40
    local col4 = 64
    
    // Row 1: Acceptable violation threshold M and violation measurement mode
    local M_fmt : di %9.3f `threshold'
    local M_fmt = strtrim("`M_fmt'")
    if "`mode'" == "iterative" {
        local mode_display "Iterative"
    }
    else {
        local mode_display "Overall"
    }
    di as text _col(`col1') "Threshold (M):" _col(`col2') as result "`M_fmt'" ///
       as text _col(`col3') "Mode:" _col(`col4') as result "`mode_display'"
    
    // Row 2: Significance level alpha and corresponding confidence level
    local alpha_fmt : di %5.3f `alpha'
    local level_fmt : di %5.1f `level'
    di as text _col(`col1') "Alpha (α):" _col(`col2') as result "`alpha_fmt'" ///
       as text _col(`col3') "Confidence:" _col(`col4') as result "`level_fmt'%"
    
    // Row 3: p-norm exponent for severity aggregation (p=inf uses sup-norm)
    if `p' >= 1e10 {
        local p_display "∞"
    }
    else {
        local p_display : di %3.1f `p'
    }
    di as text _col(`col1') "Norm (p):" _col(`col2') as result "`p_display'" ///
       as text _col(`col3') "Simulations:" _col(`col4') as result %9.0fc `sims'
    
    // =========================================================================
    // Section 3: Time Structure
    // Display panel dimensions: total periods, treatment time, pre/post counts
    // =========================================================================
    di as text ""
    di as text "{bf:Time Structure}"
    di as text "{hline 72}"
    
    di as text _col(`col1') "Total periods (T):" _col(`col2') as result "`t'" ///
       as text _col(`col3') "Treat time (t₀):" _col(`col4') as result "`t0'"
    
    di as text _col(`col1') "Pre-treatment (T_pre):" _col(`col2') as result "`t_pre'" ///
       as text _col(`col3') "Post-treatment (T_post):" _col(`col4') as result "`t_post'"
    
    local n_fmt : di %9.0fc `n'
    local n_fmt = strtrim("`n_fmt'")
    di as text _col(`col1') "Sample size (n):" _col(`col2') as result "`n_fmt'"
    
    // =========================================================================
    // Section 4: Pre-Test Results
    // Display test statistics: estimated severity, kappa constant, critical value
    // The pre-test determines whether extrapolation is justified (phi=0 => PASS)
    // =========================================================================
    di as text ""
    di as text "{bf:Pre-Test Results}"
    di as text "{hline 72}"
    
    // Estimated severity of pre-treatment parallel trends violations
    // S_pre = (1/(T_pre-1) * sum |nu_t|^p)^{1/p} for iterative violations
    // S_pre = (1/(T_pre-1) * sum |nu_bar_t|^p)^{1/p} for overall violations
    local S_pre_fmt : di %9.4f `s_pre'
    local S_pre_fmt = strtrim("`S_pre_fmt'")
    di as text _col(`col1') "Severity (Ŝ_pre):" _col(`col2') as result "`S_pre_fmt'"
    
    // Kappa: bias amplification factor from Proposition 1 (Mikhaeil & Harshaw 2025)
    // kappa = (1/T_post * sum_{t=1}^{T_post} t^q)^{1/q} where 1/p + 1/q = 1
    // Only applicable in iterative mode; overall mode bounds bias directly by S_pre
    if "`mode'" == "iterative" {
        local kappa_fmt : di %9.3f `kappa'
        local kappa_fmt = strtrim("`kappa_fmt'")
        di as text _col(`col1') "Kappa (κ):" _col(`col2') as result "`kappa_fmt'"
    }
    else {
        di as text _col(`col1') "Kappa (κ):" _col(`col2') as result "{it:N/A}" as text "  (overall mode uses Ŝ_pre directly)"
    }
    
    // Critical value from the (1-alpha) quantile of psi(Z) where Z ~ N(0, Sigma_hat)
    // See Theorem 2 and Appendix D.5 in Mikhaeil & Harshaw (2025)
    local f_alpha_fmt : di %9.3f `f_alpha'
    local f_alpha_fmt = strtrim("`f_alpha_fmt'")
    di as text _col(`col1') "Critical value f(α,Σ̂):" _col(`col2') as result "`f_alpha_fmt'"
    
    // Pre-test decision based on Theorem 1: phi = 1{S_pre_hat > M}
    // phi = 0 => PASS: extrapolation condition satisfied, proceed with inference
    // phi = 1 => FAIL: extrapolation condition violated, inference not justified
    // phi = . => Data issue preventing computation
    di as text ""
    if missing(`phi') {
        di as text _col(`col1') "Pre-test result:" ///
           _col(`col2') as error "{bf:{col 32}N/A}" ///
           as text " (Data issue - missing observations in some periods)"
    }
    else if `phi' == 0 {
        di as text _col(`col1') "Pre-test result:" ///
           _col(`col2') as result "{bf:{col 32}PASS}" ///
           as text " (Ŝ_pre ≤ M, extrapolation acceptable)"
    }
    else {
        di as text _col(`col1') "Pre-test result:" ///
           _col(`col2') as text "{bf:{col 32}FAIL}" ///
           as text " (Ŝ_pre > M, extrapolation rejected)"
    }
    
    // =========================================================================
    // Section 5: Treatment Effect Estimates
    // Display DID estimate and confidence intervals (conditional and conventional)
    // =========================================================================
    di as text ""
    di as text "{bf:Treatment Effect Estimates}"
    di as text "{hline 72}"
    
    // Average post-treatment DID estimate: delta_bar = (1/T_post) * sum delta_t
    // This estimates the identified DID estimand, which equals ATT under parallel trends
    _pretest_format_number `delta_bar'
    local delta_bar_fmt = r(formatted)
    di as text _col(`col1') "δ̄ (change vs t₀):" _col(`col2') as result "`delta_bar_fmt'"
    
    // Conditionally valid confidence interval from Theorem 2
    // Valid only when phi = 0 (pre-test passes) under well-separated null
    if missing(`phi') {
        di as text _col(`col1') "Conditional CI:" ///
           _col(`col2') as result "{it:[Not available - data issue]}"
    }
    else if `phi' == 0 {
        // Pre-test passed: display conditionally valid CI for ATT
        if `ci_lower' != -999 & `ci_upper' != -999 {
            _pretest_format_number `ci_lower'
            local ci_lo_fmt = r(formatted)
            _pretest_format_number `ci_upper'
            local ci_hi_fmt = r(formatted)
            di as text _col(`col1') "Conditional CI (`level'%):" ///
               _col(`col2') as result "[`ci_lo_fmt', `ci_hi_fmt']"
        }
        
        // Conventional CI for comparison (assumes exact parallel trends)
        if `ci_conv_lower' != -999 & `ci_conv_upper' != -999 & !missing(`ci_conv_lower') & !missing(`ci_conv_upper') {
            _pretest_format_number `ci_conv_lower'
            local conv_lo_fmt = r(formatted)
            _pretest_format_number `ci_conv_upper'
            local conv_hi_fmt = r(formatted)
            di as text _col(`col1') "Conventional CI (`level'%):" ///
               _col(`col2') as result "[`conv_lo_fmt', `conv_hi_fmt']" ///
               as text "  {it:(assuming parallel trends)}"
        }
    }
    else {
        // Pre-test failed: conditionally valid CI is not justified
        di as text _col(`col1') "Conditional CI:" ///
           _col(`col2') as result "{it:[Not available - pre-test failed]}"
        
        // Conventional CI shown for reference but interpretation requires caution
        if `ci_conv_lower' != -999 & `ci_conv_upper' != -999 & !missing(`ci_conv_lower') & !missing(`ci_conv_upper') {
            _pretest_format_number `ci_conv_lower'
            local conv_lo_fmt = r(formatted)
            _pretest_format_number `ci_conv_upper'
            local conv_hi_fmt = r(formatted)
            di as text _col(`col1') "Conventional CI (`level'%):" ///
               _col(`col2') as result "[`conv_lo_fmt', `conv_hi_fmt']" ///
               as text "  {it:(assuming parallel trends)}"
        }
    }
    
    // =========================================================================
    // Section 6: Footer
    // Methodological notes explaining the test statistic and CI formulas
    // =========================================================================
    di as text ""
    di as text "{hline 72}"
    di as text "Note: φ = 1{Ŝ_pre > M}; PASS means φ=0 (extrapolation justified)"
    if "`mode'" == "iterative" {
        di as text "      Iterative mode: CI = δ̄ ± {κ·Ŝ_pre + f(α,Σ̂)/√n}"
    }
    else {
        di as text "      Overall mode: CI = δ̄ ± {Ŝ_pre + f(α,Σ̂)/√n} (no κ)"
    }
    di as text "{hline 72}"
    di as text "{it:Note: δ̄ measures change relative to t₀, NOT the ATT level. See {help pretest##results:help pretest}.}"
end

// =============================================================================
// Helper Program: Numeric Formatting for Display
// =============================================================================

*! _pretest_format_number
*! 
*! Title: Adaptive Numeric Formatter for Econometric Output
*!
*! Description:
*!   Selects appropriate display format based on value magnitude to ensure
*!   readable output across the range of values encountered in DID analysis.
*!
*! Formatting Rules:
*!   Missing values:        "."
*!   Zero:                  "0.000"
*!   Small (|x| < 0.0005):  %9.6f (six decimal places)
*!   Standard:              %9.3f (three decimal places)
*!   Large (|x| >= 1000):   %12.0fc (comma-separated)
*!   Extreme (|x| >= 1e6):  %9.3e (scientific notation)
*!
*! Returns:
*!   r(formatted) - formatted string representation of the input value

program define _pretest_format_number, rclass
    version 17.0
    args value
    
    if missing(`value') {
        return local formatted "."
        exit
    }
    
    local abs_val = abs(`value')
    
    if `abs_val' == 0 {
        return local formatted "0.000"
        exit
    }
    
    // Scientific notation for very large values
    if `abs_val' >= 1e6 {
        local fmt : di %9.3e `value'
        local fmt = strtrim("`fmt'")
        return local formatted "`fmt'"
        exit
    }
    
    // Increased precision for small values that would round to zero
    if `abs_val' > 0 & `abs_val' < 0.0005 {
        local fmt : di %9.6f `value'
        local fmt = strtrim("`fmt'")
        return local formatted "`fmt'"
        exit
    }
    
    // Comma-separated format for large values
    if `abs_val' >= 1000 {
        local fmt : di %12.0fc `value'
        local fmt = strtrim("`fmt'")
        return local formatted "`fmt'"
        exit
    }
    
    // Standard three-decimal format
    local fmt : di %9.3f `value'
    local fmt = strtrim("`fmt'")
    return local formatted "`fmt'"
end
