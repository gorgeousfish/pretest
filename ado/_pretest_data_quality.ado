*! _pretest_data_quality.ado v0.1.0
*! Data Quality Validation for Difference-in-Differences Pre-Test Analysis
*!
*! Purpose:
*!   Validates that panel data satisfy the structural requirements for
*!   implementing the conditional extrapolation pre-test methodology
*!   proposed by Mikhaeil and Harshaw (2025).
*!
*! The validation ensures:
*!   - Block-adoption design: binary treatment indicator D in {0,1}
*!   - Sufficient time periods: T >= 3 (required for T_pre >= 2)
*!   - Balanced panel structure for consistent mean estimation
*!   - Adequate variation in treatment and outcome variables
*!
*! Syntax:
*!   _pretest_data_quality outcome, treatment(varname) time(varname) [options]
*!
*! Required:
*!   outcome            - Outcome variable Y_t
*!   treatment(varname) - Binary treatment group indicator D in {0,1}
*!   time(varname)      - Time period variable t = 1, ..., T
*!
*! Options:
*!   cluster(varname)   - Clustering variable for robust inference
*!   strict             - Fail if any quality issue is detected
*!   noreport           - Suppress detailed diagnostic output
*!
*! Returns:
*!   r(quality_passed)  - 1 if all critical checks passed, 0 otherwise
*!   r(n_issues)        - Count of detected quality issues
*!   r(T)               - Number of distinct time periods
*!   r(balanced_panel)  - 1 if panel is balanced, 0 otherwise, . if unknown
*!   r(has_variation)   - 1 if treatment has cross-sectional variation
*!   r(treat_rate)      - Proportion of observations in treatment group
*!   r(n_missing)       - Count of missing values in key variables
*!   r(n_clusters)      - Number of clusters (if cluster option specified)
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test:
*!   Valid Inference when Testing Violations of Parallel Trends for
*!   Difference-in-Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

program define _pretest_data_quality, rclass
    version 17.0
    
    syntax varlist(min=1 max=1 numeric) , ///
        TREATment(varname numeric) ///
        TIME(varname numeric) ///
        [CLuster(varname) ///
         STRICT ///
         NOReport]
    
    local outcome `varlist'
    
    local n_issues = 0
    local quality_passed = 1
    local treat_rate = .
    local is_balanced = .
    
    if "`noreport'" == "" {
        di as text _n "{hline 70}"
        di as text "Data Quality Validation for Pre-Test Analysis"
        di as text "{hline 70}"
    }
    
    // =========================================================================
    // CHECK 1: Missing Values
    // Rationale: Complete cases required for consistent estimation of group
    //            means E[Y_t | D=d] used in constructing delta_t and nu_t.
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "CHECK 1: Missing Values"
    }
    
    quietly count if missing(`outcome')
    local n_miss_outcome = r(N)
    quietly count if missing(`treatment')
    local n_miss_treat = r(N)
    quietly count if missing(`time')
    local n_miss_time = r(N)
    
    local total_missing = `n_miss_outcome' + `n_miss_treat' + `n_miss_time'
    
    if `total_missing' > 0 {
        local n_issues = `n_issues' + 1
        if "`noreport'" == "" {
            di as text "  WARNING: Missing values detected"
            if `n_miss_outcome' > 0 {
                di as text "    `outcome' (outcome): " as result %8.0f `n_miss_outcome' as text " missing"
            }
            if `n_miss_treat' > 0 {
                di as text "    `treatment' (treatment): " as result %8.0f `n_miss_treat' as text " missing"
            }
            if `n_miss_time' > 0 {
                di as text "    `time' (time): " as result %8.0f `n_miss_time' as text " missing"
            }
        }
        if "`strict'" != "" {
            local quality_passed = 0
        }
    }
    else if "`noreport'" == "" {
        di as text "  PASSED: No missing values in key variables"
    }
    
    // =========================================================================
    // CHECK 2: Binary Treatment Indicator
    // Rationale: The pre-test methodology assumes block-adoption design where
    //            D in {0,1} indicates treatment group membership. This ensures
    //            proper construction of group means for delta_t and nu_t.
    //            See Mikhaeil & Harshaw (2025), Section 2.1.
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "CHECK 2: Treatment Variable Binary"
    }
    
    quietly levelsof `treatment', local(treat_levels)
    local n_treat_levels : word count `treat_levels'
    
    local is_binary = 1
    foreach val of local treat_levels {
        if `val' != 0 & `val' != 1 {
            local is_binary = 0
        }
    }
    
    if `n_treat_levels' != 2 | `is_binary' == 0 {
        local n_issues = `n_issues' + 1
        if "`noreport'" == "" {
            di as text "  WARNING: Treatment not binary (0/1)"
            di as text "    Unique values: `treat_levels'"
        }
        if "`strict'" != "" {
            local quality_passed = 0
        }
    }
    else if "`noreport'" == "" {
        di as text "  PASSED: Treatment is binary (0/1)"
    }
    
    // =========================================================================
    // CHECK 3: Treatment Group Variation
    // Rationale: Both treatment (D=1) and control (D=0) groups must have
    //            observations to compute the DID estimand:
    //            delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
    //            See Mikhaeil & Harshaw (2025), Section 2.1.
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "CHECK 3: Treatment Variation"
    }
    
    quietly summarize `treatment'
    local treat_mean = r(mean)
    local has_variation = (`treat_mean' > 0 & `treat_mean' < 1)
    
    if `has_variation' == 0 {
        local n_issues = `n_issues' + 1
        if "`noreport'" == "" {
            di as text "  ERROR: No treatment variation (all treated or all control)"
            di as text "    Treatment mean: " as result %9.4f `treat_mean'
        }
        local quality_passed = 0
    }
    else {
        quietly count if `treatment' == 1
        local n_treated = r(N)
        quietly count
        local n_total = r(N)
        local treat_rate = `n_treated' / `n_total'
        
        if "`noreport'" == "" {
            di as text "  PASSED: Treatment has variation"
            di as text "    Treatment rate: " as result %9.4f `treat_rate'
            di as text "    N treated: " as result %8.0f `n_treated' as text " / " as result %8.0f `n_total'
        }
    }
    
    // =========================================================================
    // CHECK 4: Sufficient Time Periods (T >= 3)
    // Rationale: The pre-test requires T_pre >= 2 to estimate at least one
    //            iterative violation nu_t = E[Y_t - Y_{t-1} | D=1] -
    //            E[Y_t - Y_{t-1} | D=0] for t in {2, ..., t0-1}.
    //            With at least one post-treatment period, this requires T >= 3.
    //            See Mikhaeil & Harshaw (2025), Section 2.1.
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "CHECK 4: Time Periods (T >= 3 required)"
    }
    
    quietly levelsof `time', local(time_levels)
    local T : word count `time_levels'
    
    if `T' < 3 {
        local n_issues = `n_issues' + 1
        if "`noreport'" == "" {
            di as text "  ERROR: Insufficient time periods"
            di as text "    T = " as result %3.0f `T' as text " (minimum required: 3)"
        }
        local quality_passed = 0
    }
    else if "`noreport'" == "" {
        di as text "  PASSED: Sufficient time periods"
        di as text "    T = " as result %3.0f `T'
    }
    
    // =========================================================================
    // CHECK 5: Panel Balance
    // Rationale: Balanced panels ensure consistent sample composition across
    //            time periods when computing group means. Unbalanced panels
    //            may introduce composition bias in the estimators.
    //            The methodology assumes n_t observations per period.
    //            See Mikhaeil & Harshaw (2025), Section 2.3.
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "CHECK 5: Panel Balance"
    }
    
    capture xtset
    if _rc == 0 {
        local panelvar = r(panelvar)
        local timevar = r(timevar)
        
        tempvar panel_count
        quietly egen `panel_count' = count(`outcome'), by(`panelvar')
        quietly summarize `panel_count'
        local min_T = r(min)
        local max_T = r(max)
        local is_balanced = (`min_T' == `max_T')
        
        if `is_balanced' == 0 {
            local n_issues = `n_issues' + 1
            if "`noreport'" == "" {
                di as text "  WARNING: Unbalanced panel"
                di as text "    Min obs per unit: " as result %4.0f `min_T'
                di as text "    Max obs per unit: " as result %4.0f `max_T'
            }
        }
        else if "`noreport'" == "" {
            di as text "  PASSED: Balanced panel"
            di as text "    Obs per unit: " as result %4.0f `min_T'
        }
    }
    else {
        local is_balanced = .
        if "`noreport'" == "" {
            di as text "  SKIPPED: Panel not set (use xtset first)"
        }
    }
    
    // =========================================================================
    // CHECK 6: Outcome Variable Variation
    // Rationale: Outcome variation is necessary for meaningful estimation of
    //            treatment effects. A constant outcome implies tau_t = 0 for
    //            all t, making the pre-test analysis uninformative.
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "CHECK 6: Outcome Variable Variation"
    }
    
    quietly summarize `outcome'
    local outcome_sd = r(sd)
    local outcome_min = r(min)
    local outcome_max = r(max)
    
    if `outcome_sd' == 0 | `outcome_sd' == . {
        local n_issues = `n_issues' + 1
        if "`noreport'" == "" {
            di as text "  ERROR: No variation in outcome variable"
        }
        local quality_passed = 0
    }
    else if "`noreport'" == "" {
        di as text "  PASSED: Outcome has variation"
        di as text "    Range: [" as result %12.4f `outcome_min' as text ", " as result %12.4f `outcome_max' as text "]"
        di as text "    Std. Dev.: " as result %12.4f `outcome_sd'
    }
    
    // =========================================================================
    // CHECK 7: Cluster Variable (if specified)
    // Rationale: Cluster-robust inference requires a sufficient number of
    //            clusters for asymptotic validity. Standard recommendations
    //            suggest G >= 30-50 clusters, though this is context-dependent.
    //            With few clusters, inference may be unreliable.
    //            See Cameron & Miller (2015) for cluster-robust inference.
    // =========================================================================
    if "`cluster'" != "" {
        if "`noreport'" == "" {
            di as text _n "CHECK 7: Cluster Variable"
        }
        
        quietly levelsof `cluster', local(cluster_levels)
        local n_clusters : word count `cluster_levels'
        
        if `n_clusters' < 10 {
            local n_issues = `n_issues' + 1
            if "`noreport'" == "" {
                di as text "  WARNING: Few clusters (may affect inference)"
                di as text "    Number of clusters: " as result %6.0f `n_clusters'
            }
        }
        else if "`noreport'" == "" {
            di as text "  PASSED: Sufficient clusters"
            di as text "    Number of clusters: " as result %6.0f `n_clusters'
        }
        
        return scalar n_clusters = `n_clusters'
    }
    
    // =========================================================================
    // SUMMARY
    // =========================================================================
    if "`noreport'" == "" {
        di as text _n "{hline 70}"
        di as text "SUMMARY"
        di as text "{hline 70}"
        
        if `quality_passed' == 1 & `n_issues' == 0 {
            di as result "  All quality checks passed."
            di as text "  Data structure is suitable for pre-test analysis."
        }
        else if `quality_passed' == 1 {
            di as text "  Quality check passed with " as result %2.0f `n_issues' as text " warning(s)."
            di as text "  Review warnings before proceeding with analysis."
        }
        else {
            di as error "  Quality check FAILED with " as result %2.0f `n_issues' as text " issue(s)."
            di as text "  Address critical issues before proceeding."
        }
        
        di as text "{hline 70}"
    }
    
    // =========================================================================
    // Return Results
    // =========================================================================
    return scalar quality_passed = `quality_passed'
    return scalar n_issues = `n_issues'
    return scalar T = `T'
    return scalar balanced_panel = `is_balanced'
    return scalar has_variation = `has_variation'
    return scalar treat_rate = `treat_rate'
    return scalar n_missing = `total_missing'
    
end
