*! _pretest_output v0.1.0 - Formatted Output Display
*! 
*! Description:
*!   Displays formatted results for the pretest command following
*!   standard Stata output conventions.
*!
*! Output sections:
*!   1. Header     - Command name, reference
*!   2. Parameters - Threshold M, mode, significance level
*!   3. Time       - T, t0, T_pre, T_post
*!   4. Pre-test   - PASS/FAIL, severity, kappa
*!   5. Estimates  - ATT and confidence intervals
*!   6. Footer     - Notes and formulas
*!
*! Numeric formatting:
*!   - Standard values: %9.3f
*!   - Large numbers (>=1000): %9.0fc (comma-separated)
*!   - Extreme values (|x|>=1e6 or |x|<1e-4): %9.3e (scientific)
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), "In Defense of the Pre-Test"
*!   arXiv:2510.26470

program define _pretest_output, rclass
    version 17.0
    
    // Note: Uses rclass to avoid clearing caller's e() results
    
    // PHI is real type to support missing values (.)
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
    
    // ========================================
    // Section 1: Header
    // ========================================
    di as text ""
    di as text "{hline 72}"
    di as text "{bf:Conditional Extrapolation Pre-Test for Difference-in-Differences}"
    di as text "Mikhaeil & Harshaw (2025)"
    di as text "{hline 72}"
    
    // ========================================
    // Section 2: Model Specification
    // ========================================
    di as text ""
    di as text "{bf:Model Specification}"
    di as text "{hline 72}"
    
    // Column alignment for parameter table
    local col1 = 4
    local col2 = 28
    local col3 = 40
    local col4 = 64
    
    // Row 1: Threshold M and Mode
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
    
    // Row 2: Significance level and confidence level
    local alpha_fmt : di %5.3f `alpha'
    local level_fmt : di %5.1f `level'
    di as text _col(`col1') "Alpha (α):" _col(`col2') as result "`alpha_fmt'" ///
       as text _col(`col3') "Confidence:" _col(`col4') as result "`level_fmt'%"
    
    // Row 3: Norm parameter p
    if `p' >= 1e10 {
        local p_display "∞"
    }
    else {
        local p_display : di %3.1f `p'
    }
    di as text _col(`col1') "Norm (p):" _col(`col2') as result "`p_display'" ///
       as text _col(`col3') "Simulations:" _col(`col4') as result %9.0fc `sims'
    
    // ========================================
    // Section 3: Time Structure
    // ========================================
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
    
    // ========================================
    // Section 4: Pre-Test Results
    // ========================================
    di as text ""
    di as text "{bf:Pre-Test Results}"
    di as text "{hline 72}"
    
    // Severity estimate (inline formatting)
    local S_pre_fmt : di %9.4f `s_pre'
    local S_pre_fmt = strtrim("`S_pre_fmt'")
    di as text _col(`col1') "Severity (Ŝ_pre):" _col(`col2') as result "`S_pre_fmt'"
    
    // Kappa constant (iterative mode only) or N/A for overall mode
    if "`mode'" == "iterative" {
        local kappa_fmt : di %9.3f `kappa'
        local kappa_fmt = strtrim("`kappa_fmt'")
        di as text _col(`col1') "Kappa (κ):" _col(`col2') as result "`kappa_fmt'"
    }
    else {
        di as text _col(`col1') "Kappa (κ):" _col(`col2') as result "{it:N/A}" as text "  (overall mode uses Ŝ_pre directly)"
    }
    
    // Critical value f(α, Σ̂)
    local f_alpha_fmt : di %9.3f `f_alpha'
    local f_alpha_fmt = strtrim("`f_alpha_fmt'")
    di as text _col(`col1') "Critical value f(α,Σ̂):" _col(`col2') as result "`f_alpha_fmt'"
    
    // PASS/FAIL result
    // φ = 1{Ŝ_pre > M}: φ=0 means PASS, φ=1 means FAIL, φ=. indicates data issue
    di as text ""
    if missing(`phi') {
        // Data issue: cannot compute
        di as text _col(`col1') "Pre-test result:" ///
           _col(`col2') as error "{bf:{col 32}N/A}" ///
           as text " (Data issue - missing observations in some periods)"
    }
    else if `phi' == 0 {
        // PASS
        di as text _col(`col1') "Pre-test result:" ///
           _col(`col2') as result "{bf:{col 32}PASS}" ///
           as text " (Ŝ_pre ≤ M, extrapolation acceptable)"
    }
    else {
        // FAIL
        di as text _col(`col1') "Pre-test result:" ///
           _col(`col2') as error "{bf:{col 32}FAIL}" ///
           as text " (Ŝ_pre > M, extrapolation rejected)"
    }
    
    // ========================================
    // Section 5: Treatment Effect Estimates
    // ========================================
    di as text ""
    di as text "{bf:Treatment Effect Estimates}"
    di as text "{hline 72}"
    
    // ATT estimate (inline formatting)
    // Note: δ̄ measures change relative to t₀
    _pretest_format_number `delta_bar'
    local delta_bar_fmt = r(formatted)
    di as text _col(`col1') "δ̄ (change vs t₀):" _col(`col2') as result "`delta_bar_fmt'"
    
    // Conditionally valid CI (only shown when pretest passes)
    if missing(`phi') {
        // Data issue: CI not available
        di as text _col(`col1') "Conditional CI:" ///
           _col(`col2') as result "{it:[Not available - data issue]}"
    }
    else if `phi' == 0 {
        // PASS: show conditionally valid CI
        if `ci_lower' != -999 & `ci_upper' != -999 {
            _pretest_format_number `ci_lower'
            local ci_lo_fmt = r(formatted)
            _pretest_format_number `ci_upper'
            local ci_hi_fmt = r(formatted)
            di as text _col(`col1') "Conditional CI (`level'%):" ///
               _col(`col2') as result "[`ci_lo_fmt', `ci_hi_fmt']"
        }
        
        // Conventional CI (if provided and non-missing)
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
        // FAIL: show CI not available message
        di as text _col(`col1') "Conditional CI:" ///
           _col(`col2') as result "{it:[Not available - pre-test failed]}"
        
        // Conventional CI still shown (if provided and non-missing)
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
    
    // ========================================
    // Section 6: Footer
    // ========================================
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

// ============================================================================
// Helper Function: Numeric Formatting
// ============================================================================

/*
 * _pretest_format_number - Select appropriate format based on value magnitude
 *
 * Formatting rules:
 *   1. Standard values: %9.3f
 *   2. Large numbers (|x| >= 1000): %9.0fc (comma-separated)
 *   3. Extreme values (|x| >= 1e6 or 0 < |x| < 1e-4): %9.3e (scientific)
 *   4. Zero: "0.000"
 *   5. Missing: "."
 */
program define _pretest_format_number, rclass
    version 17.0
    args value
    
    // Missing value
    if missing(`value') {
        return local formatted "."
        exit
    }
    
    local abs_val = abs(`value')
    
    // Zero value
    if `abs_val' == 0 {
        return local formatted "0.000"
        exit
    }
    
    // Extreme large value (|x| >= 1e6): scientific notation
    if `abs_val' >= 1e6 {
        local fmt : di %9.3e `value'
        local fmt = strtrim("`fmt'")
        return local formatted "`fmt'"
        exit
    }
    
    // Small value that would display as 0.000 or -0.000 (0 < |x| < 0.0005): use more precision
    // Threshold 0.0005 ensures values that round to 0.000 at 3 decimals get better display
    if `abs_val' > 0 & `abs_val' < 0.0005 {
        // Use 6 decimal places for small values
        local fmt : di %9.6f `value'
        local fmt = strtrim("`fmt'")
        return local formatted "`fmt'"
        exit
    }
    
    // Large number (|x| >= 1000): comma-separated
    if `abs_val' >= 1000 {
        local fmt : di %12.0fc `value'
        local fmt = strtrim("`fmt'")
        return local formatted "`fmt'"
        exit
    }
    
    // Standard value: %9.3f
    local fmt : di %9.3f `value'
    local fmt = strtrim("`fmt'")
    return local formatted "`fmt'"
end

