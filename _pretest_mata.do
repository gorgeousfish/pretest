*! _pretest_mata.do
*! Mata Library Loader for the pretest Package
*! Version: 0.1.0
*!
*! Purpose:
*!   Loads all Mata functions required by the pretest command.
*!   This script is automatically invoked by pretest.ado or can be
*!   executed manually for development and testing purposes.
*!
*! Usage:
*!   do _pretest_mata.do
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test:
*!   Valid Inference when Testing Violations of Parallel Trends for
*!   Difference-in-Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

// Display initialization status
di as text "{hline 60}"
di as text "Initializing pretest Mata library..."
di as text "{hline 60}"

// Clear existing Mata function definitions to prevent version conflicts
capture noisily mata: mata drop _pretest_*()

// ============================================================================
// Mata Source File Loader
// ----------------------------------------------------------------------------
// Loads Mata source files in dependency order using multiple search strategies
// for compatibility with Stata's package installation conventions.
// ============================================================================

// Helper program: Locate and execute Mata source files
// Search priority:
//   1. mata/ subdirectory relative to c(filename)
//   2. mata/ subdirectory relative to findfile result for this script
//   3. Direct findfile search for the target file
//   4. mata/ subdirectory relative to pretest.ado location
//   5. Standard sysdir locations (PLUS, PERSONAL, SITE)
capture program drop _pretest_load_mata
program define _pretest_load_mata
    args filename
    
    // Strategy 1: Locate via current script directory (c(filename))
    local thisfile "`c(filename)'"
    if "`thisfile'" != "" {
        local basedir = subinstr("`thisfile'", "_pretest_mata.do", "", .)
        local matapath "`basedir'mata/`filename'"
        capture confirm file "`matapath'"
        if _rc == 0 {
            quietly do "`matapath'"
            di as text "  [OK] `filename' loaded from mata/"
            exit
        }
    }
    
    // Strategy 2: Locate via findfile on this loader script
    capture findfile _pretest_mata.do
    if _rc == 0 {
        local basedir = subinstr("`r(fn)'", "_pretest_mata.do", "", .)
        local matapath "`basedir'mata/`filename'"
        capture confirm file "`matapath'"
        if _rc == 0 {
            quietly do "`matapath'"
            di as text "  [OK] `filename' loaded from mata/"
            exit
        }
    }
    
    // Strategy 3: Direct findfile search for target file
    capture findfile `filename'
    if _rc == 0 {
        quietly do "`r(fn)'"
        di as text "  [OK] `filename' loaded via findfile"
        exit
    }
    
    // Strategy 4: Locate via pretest.ado sibling mata/ directory
    capture findfile pretest.ado
    if _rc == 0 {
        local adomatch "`r(fn)'"
        local basedir = subinstr("`adomatch'", "pretest.ado", "", .)
        local matapath "`basedir'mata/`filename'"
        capture confirm file "`matapath'"
        if _rc == 0 {
            quietly do "`matapath'"
            di as text "  [OK] `filename' loaded via pretest.ado sibling"
            exit
        }
    }
    
    // Strategy 5: Search standard Stata system directories
    foreach sdir in PLUS PERSONAL SITE {
        local sysdir : sysdir `sdir'
        if "`sysdir'" != "" {
            // Try direct in sysdir/mata/
            local matapath "`sysdir'mata/`filename'"
            capture confirm file "`matapath'"
            if _rc == 0 {
                quietly do "`matapath'"
                di as text "  [OK] `filename' loaded from `sdir'/mata/"
                exit
            }
        }
    }
    
    di as error "  [X] `filename' not found"
    di as error "      Searched: adopath, c(filename), findfile, pretest.ado sibling, sysdir"
    exit 601
end

// ----------------------------------------------------------------------------
// Load Mata modules in dependency order
// ----------------------------------------------------------------------------

// 1. Utility functions: Matrix operations, numerical helpers (no dependencies)
_pretest_load_mata _pretest_utils.mata

// 2. DID estimators: delta_t computation per Section 2.1
_pretest_load_mata _pretest_estimators.mata

// 3. Violation estimators: nu_t and S_pre per Section 3.1
_pretest_load_mata _pretest_violations.mata

// 4. Covariance estimation: Sigma matrix per Assumption 2
_pretest_load_mata _pretest_covariance.mata

// 5. Kappa computation: Bias bound constant per Proposition 1
_pretest_load_mata _pretest_kappa.mata

// 6. Psi function: Critical value computation per Section 5.1
_pretest_load_mata _pretest_psi.mata

// 7. Confidence intervals: Implementation of Theorem 2
_pretest_load_mata _pretest_ci.mata

// 8. Simulation module: Coverage analysis (optional, for development)
capture _pretest_load_mata _pretest_simulation.mata
if _rc {
    di as text "  [!] _pretest_simulation.mata not loaded (optional module)"
}

// 9. Main computation: Orchestrates full pretest analysis
_pretest_load_mata _pretest_main.mata

// Clean up temporary helper program
capture program drop _pretest_load_mata

// ============================================================================
// Initialization Complete
// ============================================================================

di as text "{hline 60}"
di as text "pretest Mata library initialized successfully."
di as text "{hline 60}"
