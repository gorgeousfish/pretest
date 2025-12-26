*! pretest_coverage.ado
*! Coverage Rate Simulation for Pretest Framework
*! Monte Carlo simulation to validate conditional coverage
*!
*! Syntax:
*!   pretest_coverage, [options]
*!
*! Options:
*!   nsims(#)      - Number of simulation replications (default: 500)
*!   nobs(#)       - Number of units per simulation (default: 200)
*!   periods(#)    - Number of time periods T (default: 5)
*!   t0(#)         - Treatment time (default: 3)
*!   tau(#)        - True ATT effect (default: 1)
*!   ptreat(#)     - Treatment probability (default: 0.3)
*!   threshold(#)  - Pretest threshold M (default: 0.5)
*!   p(#)          - Severity p-norm (default: 2)
*!   level(#)      - Confidence level (default: 95)
*!   mode(string)  - "iterative" or "overall" (default: iterative)
*!   violation(string) - Violation pattern: "zero", "constant", "linear" (default: zero)
*!   vscale(#)     - Violation scale factor (default: 0.1)
*!   mcsims(#)     - MC sims for critical value per replication (default: 1000)
*!   seed(#)       - Base random seed (default: 12345)
*!
*! Returns (r-class):
*!   r(coverage_rate)   - Conditional coverage rate
*!   r(pass_rate)       - Pretest pass rate
*!   r(effective_coverage) - Pass rate × Coverage
*!   r(n_passed)        - Number of passed simulations
*!   r(n_sims)          - Total simulations
*!   r(mean_ci_width)   - Average CI width

program define pretest_coverage, rclass
    version 17.0
    
    syntax , ///
        [NSims(integer 500) ///
         NOBS(integer 200) ///
         PERiods(integer 5) ///
         T0(integer 3) ///
         TAU(real 1) ///
         PTReat(real 0.3) ///
         THReshold(real 0.5) ///
         P(real 2) ///
         Level(cilevel) ///
         MODE(string) ///
         VIOLation(string) ///
         VSCale(real 0.1) ///
         MCSims(integer 1000) ///
         SEED(integer 12345)]
    
    // Default mode
    if "`mode'" == "" {
        local mode "iterative"
    }
    
    // Validate mode
    if "`mode'" != "iterative" & "`mode'" != "overall" {
        di as error "Error: mode must be 'iterative' or 'overall'"
        exit 198
    }
    
    // Default violation pattern
    if "`violation'" == "" {
        local violation "zero"
    }
    
    // Validate violation pattern
    if "`violation'" != "zero" & "`violation'" != "constant" & ///
       "`violation'" != "linear" & "`violation'" != "quadratic" {
        di as error "Error: violation must be 'zero', 'constant', 'linear', or 'quadratic'"
        exit 198
    }
    
    // Calculate alpha
    local alpha = 1 - `level'/100
    
    // Calculate T_pre
    local T_pre = `t0' - 1
    
    // Display header
    di as text _n "{hline 78}"
    di as text "{center:PRETEST Coverage Rate Simulation}"
    di as text "{center:Monte Carlo Validation for Mikhaeil & Harshaw (2025)}"
    di as text "{hline 78}"
    
    // Display parameters
    di as text _n "Simulation Parameters:"
    di as text "  Replications:         " as result %8.0f `nsims'
    di as text "  Units per sim (n):    " as result %8.0f `nobs'
    di as text "  Time periods (T):     " as result %8.0f `periods'
    di as text "  Treatment time (t0):  " as result %8.0f `t0'
    di as text "  Pre-treatment (T_pre):" as result %8.0f `T_pre'
    di as text "  True ATT (τ):         " as result %8.4f `tau'
    di as text "  Treatment prob:       " as result %8.4f `ptreat'
    di as text "  Threshold (M):        " as result %8.4f `threshold'
    di as text "  p-norm:               " as result %8.0f `p'
    di as text "  Mode:                 " as result "`mode'"
    di as text "  Violation pattern:    " as result "`violation'"
    di as text "  Violation scale:      " as result %8.4f `vscale'
    di as text "  MC sims per rep:      " as result %8.0f `mcsims'
    di as text "  Base seed:            " as result %8.0f `seed'
    di as text "  Confidence level:     " as result %8.1f `level' "%"
    
    // Load Mata functions
    capture findfile pretest.ado
    if _rc == 0 {
        local ado_dir = subinstr("`r(fn)'", "pretest.ado", "", .)
    }
    else {
        local ado_dir ""
    }
    quietly capture do "`ado_dir'_pretest_mata.do"
    
    // Run simulation in Mata
    di as text _n "Running simulation..."
    
    tempname cov_rate pass_rate eff_cov n_passed n_total mean_width
    
    mata: {
        real colvector nu_vec
        real rowvector results
        real scalar T_pre_m
        
        T_pre_m = `T_pre'
        
        // Generate violation vector based on pattern
        nu_vec = _pretest_violation_dgp(T_pre_m, "`violation'", `vscale')
        
        // Run coverage simulation
        results = _pretest_coverage_simulation(
            `nsims',          // n_sims
            `tau',            // tau
            `nobs',           // n
            `periods',        // T
            `t0',             // t0
            `ptreat',         // p_treat
            nu_vec,           // nu_vec
            `threshold',      // M
            `p',              // p_norm
            `alpha',          // alpha
            "`mode'",         // mode
            `mcsims',         // S_sims
            `seed'            // base_seed
        )
        
        // Store results
        st_numscalar("`cov_rate'", results[1])
        st_numscalar("`pass_rate'", results[2])
        st_numscalar("`eff_cov'", results[3])
        st_numscalar("`n_passed'", results[4])
        st_numscalar("`n_total'", results[5])
        st_numscalar("`mean_width'", results[6])
    }
    
    // Display results
    di as text _n "{hline 78}"
    di as text "{bf:Coverage Simulation Results}"
    di as text "{hline 78}"
    
    local cov = scalar(`cov_rate')
    local pass = scalar(`pass_rate')
    local eff = scalar(`eff_cov')
    local np = scalar(`n_passed')
    local nt = scalar(`n_total')
    local mw = scalar(`mean_width')
    
    di as text "  Pretest pass rate:         " as result %8.4f `pass'
    di as text "  Passed simulations:        " as result %8.0f `np' as text " / " as result %8.0f `nt'
    di as text ""
    
    // Highlight coverage rate
    if `cov' >= 0.93 {
        di as text "  Conditional coverage rate: " as result %8.4f `cov' as text " (TARGET: ≥ 0.93)"
        di as result "  ✓ Coverage target MET!"
    }
    else if !missing(`cov') {
        di as text "  Conditional coverage rate: " as text %8.4f `cov' as text " (TARGET: ≥ 0.93)"
        di as text "  ✗ Coverage target NOT met"
    }
    else {
        di as text "  Conditional coverage rate: " as result "[No simulations passed pretest]"
    }
    
    di as text ""
    di as text "  Effective coverage:        " as result %8.4f `eff'
    
    if !missing(`mw') {
        di as text "  Mean CI width:             " as result %8.4f `mw'
    }
    
    // Check theoretical targets from paper
    di as text _n "{hline 78}"
    di as text "{bf:Coverage Validation Check}"
    di as text "{hline 78}"
    
    local cov_pass = (`cov' >= 0.93 | missing(`cov'))
    local eff_pass = (`pass' >= 0.50)  // Effective reporting probability > 0.5 at minimum
    
    if `cov_pass' {
        di as text "  Conditional coverage ≥ 93%: " as result "PASS"
    }
    else {
        di as text "  Conditional coverage ≥ 93%: " as text "FAIL"
    }
    
    if `eff_pass' {
        di as text "  Effective reporting > 50%:  " as result "PASS"
    }
    else {
        di as text "  Effective reporting > 50%:  " as text "FAIL"
    }
    
    di as text "{hline 78}"
    
    // Return results
    return scalar coverage_rate = `cov'
    return scalar pass_rate = `pass'
    return scalar effective_coverage = `eff'
    return scalar n_passed = `np'
    return scalar n_sims = `nt'
    return scalar mean_ci_width = `mw'
    return scalar tau = `tau'
    return scalar threshold = `threshold'
    return scalar p = `p'
    return local mode "`mode'"
    return local violation "`violation'"
end
