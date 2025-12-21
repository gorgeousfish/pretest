*! _pretest_mata.do v0.1.0
*! Mata Library Loader for pretest Package
*!
*! Description:
*!   Loads all Mata functions required by the pretest command.
*!   Automatically called by pretest.ado or can be run manually.
*!
*! Usage:
*!   . do _pretest_mata.do
*!
*! Reference:
*!   Mikhaeil, J.M. and C. Harshaw. 2025. In Defense of the Pre-Test.
*!   arXiv:2510.26470. https://arxiv.org/abs/2510.26470

version 17.0

// Display loading status
di as text "{hline 60}"
di as text "Loading pretest Mata library..."
di as text "{hline 60}"

// Clear any existing definitions to avoid conflicts
capture noisily mata: mata drop _pretest_*()

// ============================================================
// Load Mata source files (in dependency order)
// Uses findfile for compatibility with Stata's installation structure
// ============================================================

// Helper program: locate and execute Mata files
// Priority: (0) adopath search, (1) mata/ subdirectory via c(filename), 
// (2) findfile this script, (3) findfile the target file directly
capture program drop _pretest_load_mata
program define _pretest_load_mata
    args filename
    
    // Method 0 (NEW): Search adopath for mata/ subdirectory
    // This ensures development directory takes precedence when added to adopath
    local adopaths : adopath
    local word_count : word count `adopaths'
    forvalues i = 1/`word_count' {
        local apath : word `i' of `adopaths'
        // Skip system paths
        if strpos("`apath'", "(") > 0 continue
        // Check for mata subdirectory
        local matapath "`apath'/mata/`filename'"
        capture confirm file "`matapath'"
        if _rc == 0 {
            quietly do "`matapath'"
            di as text "  [OK] `filename' loaded from adopath mata/"
            exit
        }
    }
    
    // Method 1: Locate via current script directory (c(filename))
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
    
    // Method 2: Locate via findfile on this script
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
    
    // Method 3: Direct findfile on target file
    capture findfile `filename'
    if _rc == 0 {
        quietly do "`r(fn)'"
        di as text "  [OK] `filename' loaded via findfile"
        exit
    }
    
    // Method 4: Search via findfile on pretest.ado then look in sibling mata/
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
    
    // Method 5: Search via sysdir standard locations
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

// 1. Utility functions (no dependencies)
_pretest_load_mata _pretest_utils.mata

// 2. DID estimator functions
_pretest_load_mata _pretest_estimators.mata

// 3. Violation estimation functions
_pretest_load_mata _pretest_violations.mata

// 4. Covariance matrix estimation
_pretest_load_mata _pretest_covariance.mata

// 5. Kappa constant computation
_pretest_load_mata _pretest_kappa.mata

// 6. Psi function and Monte Carlo critical values
_pretest_load_mata _pretest_psi.mata

// 7. Confidence interval functions
_pretest_load_mata _pretest_ci.mata

// 8. Coverage simulation functions (optional)
capture _pretest_load_mata _pretest_simulation.mata
if _rc {
    di as text "  [!] _pretest_simulation.mata not found (optional)"
}

// 9. Main computation function
_pretest_load_mata _pretest_main.mata

// Clean up helper program
capture program drop _pretest_load_mata

// ============================================================
// Complete
// ============================================================

di as text "{hline 60}"
di as text "pretest Mata library loaded successfully"
di as text "{hline 60}"
