*! ex_custom_plot.do - Demonstration of Custom Plotting in pretest
*! Part of the pretest-stata package
*!
*! This script demonstrates the graph customization features introduced
*! in Story 3.1, allowing researchers to tailor event study plots for
*! specific journal requirements or aesthetic preferences.
*!
*! Reference: Mikhaeil & Harshaw (2025), Figure 1

version 17.0
clear all
set more off

// ============================================================================
// SETUP: Load example data
// ============================================================================
// Load Proposition 99 dataset from package
webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"
webuse prop99_smoking.dta, clear

// Describe the data structure
describe

// ============================================================================
// EXAMPLE 1: Default Plot (Baseline - Mikhaeil-Harshaw Style)
// ============================================================================
// This produces the standard event study plot as described in the paper.
// No customization options are provided, so all defaults are used.

di _n as text "{hline 70}"
di as text "EXAMPLE 1: Default Plot (Mikhaeil-Harshaw 2025 Style)"
di as text "{hline 70}"

pretest cigsale, treatment(treated) time(year) threshold(5) treat_time(1989)

// ============================================================================
// EXAMPLE 2: Black and White Scheme for Print Journals
// ============================================================================
// Many academic journals require grayscale figures. This example demonstrates
// how to create a publication-ready black and white plot.

di _n as text "{hline 70}"
di as text "EXAMPLE 2: Black & White Scheme for Print Publications"
di as text "{hline 70}"

pretest cigsale, treatment(treated) time(year) threshold(5) treat_time(1989) ///
    ci_opt_pass(lcolor(black) lwidth(medium)) ///
    ci_opt_fail(lcolor(gs8) lpattern(dash)) ///
    line_opt_m(lcolor(gs5) lpattern(shortdash)) ///
    marker_opt_pre(mcolor(black) msymbol(O)) ///
    marker_opt_post(mcolor(black) msymbol(S)) ///
    scheme(s2mono)

// ============================================================================
// EXAMPLE 3: Suppressed Titles and Custom Labels
// ============================================================================
// For figures intended to be combined with other plots or when titles are
// added externally (e.g., in LaTeX), you may want to suppress default titles.

di _n as text "{hline 70}"
di as text "EXAMPLE 3: Suppressed Titles for LaTeX Integration"
di as text "{hline 70}"

pretest cigsale, treatment(treated) time(year) threshold(5) treat_time(1989) ///
    title("") ///
    note("") ///
    xtitle("Years Relative to Proposition 99") ///
    ytitle("Cigarette Sales (packs per capita)") ///
    legend(off)

// ============================================================================
// EXAMPLE 4: Custom Color Scheme (Journal of Finance Style)
// ============================================================================
// Different journals have different style preferences. This example shows
// a more vibrant color scheme.

di _n as text "{hline 70}"
di as text "EXAMPLE 4: Custom Color Scheme"
di as text "{hline 70}"

pretest cigsale, treatment(treated) time(year) threshold(5) treat_time(1989) ///
    ci_opt_pass(lcolor(dkgreen) lwidth(thick)) ///
    ci_opt_fail(lcolor(cranberry) lwidth(thick)) ///
    line_opt_m(lcolor(purple) lpattern(dash)) ///
    marker_opt_pre(mcolor(dknavy) msymbol(D) msize(large)) ///
    marker_opt_post(mcolor(dkgreen) msymbol(T) msize(large)) ///
    title("Custom Event Study: California Tobacco Control Program")

// ============================================================================
// EXAMPLE 5: Export to Publication Formats
// ============================================================================
// Export the graph to PNG and PDF formats for publication.

di _n as text "{hline 70}"
di as text "EXAMPLE 5: Export to PNG and PDF"
di as text "{hline 70}"

// First, generate the plot
pretest cigsale, treatment(treated) time(year) threshold(5) treat_time(1989) ///
    title("Proposition 99 Effect on Cigarette Sales")

// Export to PNG (high resolution for web/presentations)
graph export "ex_custom_plot_output.png", replace width(2400) height(1600)
di as text "Exported: ex_custom_plot_output.png"

// Export to PDF (vector format for publications)
graph export "ex_custom_plot_output.pdf", replace
di as text "Exported: ex_custom_plot_output.pdf"

// ============================================================================
// EXAMPLE 6: Overall Mode with Custom Styling
// ============================================================================
// Demonstrates that custom styling works with both iterative and overall modes.

di _n as text "{hline 70}"
di as text "EXAMPLE 6: Overall Violations Mode with Custom Style"
di as text "{hline 70}"

pretest cigsale, treatment(treated) time(year) threshold(5) treat_time(1989) ///
    overall ///
    ci_opt_pass(lcolor(teal)) ///
    marker_opt_pre(mcolor(teal) msymbol(oh)) ///
    marker_opt_post(mcolor(orange_red) msymbol(Dh))

// ============================================================================
// SUMMARY
// ============================================================================
di _n as text "{hline 70}"
di as text "Custom Plotting Examples Complete"
di as text "{hline 70}"
di as text ""
di as text "Available element-specific options:"
di as text "  ci_opt_pass()    - Confidence interval styling (pretest pass)"
di as text "  ci_opt_fail()    - Confidence interval styling (pretest fail)"
di as text "  line_opt_m()     - Threshold M line styling"
di as text "  marker_opt_pre() - Pre-treatment marker styling"
di as text "  marker_opt_post()- Post-treatment marker styling"
di as text ""
di as text "Additionally, all standard Stata twoway_options are supported:"
di as text "  title(), xtitle(), ytitle(), xlabel(), ylabel(),"
di as text "  legend(), scheme(), note(), etc."
di as text ""
di as text "For complete documentation, see: help pretest"
