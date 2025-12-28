*! build_mlib.do
*! Build pretest Mata library (.mlib) for distribution
*! 
*! Usage:
*!   cd "/path/to/pretest-stata"
*!   do build_mlib.do
*!
*! Output:
*!   lpretest.mlib - Compiled Mata library containing all pretest functions
*!
*! Note: Run this script from the pretest-stata root directory before distribution.
*!       The generated .mlib file should be included in the package.

version 17.0
clear all

// Display build information
di as text "{hline 60}"
di as text "Building pretest Mata library (lpretest.mlib)"
di as text "{hline 60}"
di as text "Working directory: `c(pwd)'"
di as text "Stata version: `c(stata_version)'"
di as text "Date: `c(current_date)' `c(current_time)'"
di as text "{hline 60}"

// Clear any existing pretest functions from memory
capture mata: mata drop _pretest_*()

// ============================================================================
// Step 1: Compile all Mata source files in dependency order
// ============================================================================

di as text _n "Step 1: Compiling Mata source files..."

// 1. Utility functions (no dependencies)
di as text "  Compiling _pretest_utils.mata..."
quietly do "mata/_pretest_utils.mata"

// 2. DID estimators
di as text "  Compiling _pretest_estimators.mata..."
quietly do "mata/_pretest_estimators.mata"

// 3. Violation estimators
di as text "  Compiling _pretest_violations.mata..."
quietly do "mata/_pretest_violations.mata"

// 4. Covariance estimation
di as text "  Compiling _pretest_covariance.mata..."
quietly do "mata/_pretest_covariance.mata"

// 5. Kappa computation
di as text "  Compiling _pretest_kappa.mata..."
quietly do "mata/_pretest_kappa.mata"

// 6. Psi function (critical value)
di as text "  Compiling _pretest_psi.mata..."
quietly do "mata/_pretest_psi.mata"

// 7. Confidence intervals
di as text "  Compiling _pretest_ci.mata..."
quietly do "mata/_pretest_ci.mata"

// 8. Simulation module (optional)
di as text "  Compiling _pretest_simulation.mata..."
capture quietly do "mata/_pretest_simulation.mata"
if _rc {
    di as text "    (skipped - optional module)"
}

// 9. Main computation engine
di as text "  Compiling _pretest_main.mata..."
quietly do "mata/_pretest_main.mata"

di as text "  All source files compiled successfully."

// ============================================================================
// Step 2: Create the Mata library (.mlib)
// ============================================================================

di as text _n "Step 2: Creating Mata library..."

// Remove existing library if present
capture erase "lpretest.mlib"

// Create new library with all _pretest_* functions
// The 'replace' option overwrites if exists
// Size 1 = small library (sufficient for this package)
mata: mata mlib create lpretest, dir(".") replace
mata: mata mlib add lpretest _pretest_*(), dir(".")

di as text "  lpretest.mlib created successfully."

// ============================================================================
// Step 3: Verify the library
// ============================================================================

di as text _n "Step 3: Verifying library contents..."

// List functions in the library
mata: mata mlib index
mata: mata describe _pretest_main_compute()

di as text _n "{hline 60}"
di as text "Build complete!"
di as text "{hline 60}"
di as text "Output file: lpretest.mlib"
di as text ""
di as text "Next steps:"
di as text "  1. Include lpretest.mlib in your package distribution"
di as text "  2. Update pretest.pkg to include the .mlib file"
di as text "  3. Users will automatically load the library via mata mlib index"
di as text "{hline 60}"

