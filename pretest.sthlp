{smcl}
{* *! version 0.1.0  22dec2025}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[XT] xtreg" "help xtreg"}{...}
{vieweralsosee "[XT] xtset" "help xtset"}{...}
{vieweralsosee "[XT] xtdidregress" "help xtdidregress"}{...}
{viewerjumpto "Syntax" "pretest##syntax"}{...}
{viewerjumpto "Description" "pretest##description"}{...}
{viewerjumpto "Options" "pretest##options"}{...}
{viewerjumpto "Examples" "pretest##examples"}{...}
{viewerjumpto "Stored results" "pretest##results"}{...}
{viewerjumpto "Methods and formulas" "pretest##methods"}{...}
{viewerjumpto "References" "pretest##references"}{...}
{viewerjumpto "Author" "pretest##author"}{...}

{title:Title}

{phang}
{bf:pretest} {hline 2} Stata module to implement the conditional extrapolation pre-test for difference-in-differences


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:pretest}
{it:outcome}
{cmd:,}
{opt treat:ment(varname)}
{opt time(varname)}
{opt thre:shold(#)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt treat:ment(varname)}}binary treatment indicator (0/1){p_end}
{synopt:{opt time(varname)}}time variable{p_end}
{synopt:{opt thre:shold(#)}}acceptable violation threshold M > 0{p_end}

{syntab:Optional}
{synopt:{opt treat_time(#)}}treatment time point t0; auto-detected only if treatment is time-varying{p_end}
{synopt:{opt p(#)}}severity norm p >= 1; default is {cmd:p(2)}{p_end}
{synopt:{opt al:pha(#)}}significance level; default is {cmd:alpha(0.05)}{p_end}
{synopt:{opt level(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt cl:uster(varname)}}cluster variable for standard errors{p_end}
{synopt:{opt over:all}}use overall violations mode (no kappa multiplier){p_end}
{synopt:{opt nog:raph}}suppress event study graph{p_end}
{synopt:{opt sim:ulate(#)}}Monte Carlo simulations; default is {cmd:simulate(5000)}{p_end}
{synopt:{opt seed(#)}}random number seed; default is {cmd:seed(12345)}{p_end}
{synopt:{opt diag:nose}}display detailed diagnostic information{p_end}

{syntab:Graph customization}
{synopt:{opt ci_opt_pass(string)}}override CI style when pretest passes{p_end}
{synopt:{opt ci_opt_fail(string)}}override CI style when pretest fails{p_end}
{synopt:{opt line_opt_m(string)}}override threshold M line style{p_end}
{synopt:{opt marker_opt_pre(string)}}override pre-treatment marker style{p_end}
{synopt:{opt marker_opt_post(string)}}override post-treatment marker style{p_end}
{synopt:{it:twoway_options}}any standard Stata {help twoway_options}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{bf:Data Requirements:}

{phang2}- {bf:Minimum 3 time periods}: T_pre >= 2 (at least two pre-treatment periods){p_end}
{phang2}- {bf:Block adoption design}: All treated units receive treatment at the same t0{p_end}
{phang2}- {bf:Binary treatment}: Treatment indicator must be 0 (control) or 1 (treated){p_end}

{pstd}
{bf:Note:} For datasets with only two time periods (e.g., pre vs post), this command
cannot be applied because iterative violations nu_t are defined only for t >= 2,
requiring at least T_pre >= 2. In such cases, use standard DID estimators:

{phang2}{cmd:. regress y i.treat##i.post}{p_end}
{phang2}{cmd:. xtdidregress (y) (d), group(id) time(t)}{p_end}

{pstd}
{cmd:pretest} implements the conditional extrapolation pre-test framework for
difference-in-differences (DID) designs proposed by Mikhaeil and Harshaw (2025).

{pstd}
Traditional pre-tests for parallel trends suffer from a fundamental problem: 
rejecting the null hypothesis of parallel trends is informative, but failing to
reject is not, since the test may simply lack power. Moreover, conditioning on 
passing such pre-tests can lead to severe bias in subsequent inference 
(Roth 2022).

{pstd}
This command provides a solution through the {it:conditional extrapolation assumption}:
if pre-treatment parallel trend violations are sufficiently small (below threshold M),
then post-treatment violations are assumed to be no worse. Under this assumption,
the command provides:

{phang2}1. A {bf:pre-test} that determines whether the estimated severity of 
pre-treatment violations falls below the user-specified threshold M.{p_end}

{phang2}2. {bf:Conditionally valid confidence intervals} for the Average Treatment 
Effect on the Treated (ATT) when the pre-test passes, with guaranteed asymptotic 
coverage conditional on passing.{p_end}

{pstd}
{bf:Key Concepts:}

{pstd}
{it:Iterative violation} (Section 2.1) at time t measures the period-to-period 
deviation from parallel trends:

{p 8 8 2}
nu_t = E[Y(0)_t - Y(0)_{t-1} | D=1] - E[Y(0)_t - Y(0)_{t-1} | D=0]

{pstd}
{it:Overall violation} (Section 2.1) at time t measures the cumulative 
deviation from parallel trends relative to the reference period t0:

{p 8 8 2}
nu_bar_t = E[Y(0)_t - Y(0)_{t0} | D=1] - E[Y(0)_t - Y(0)_{t0} | D=0]

{pstd}
Equivalently, overall violations are cumulative sums of iterative violations
(Section 2.1):

{p 8 8 2}
nu_bar_t = sum_{s=2}^{t} nu_s    where nu_bar_{t0-1} = 0 by normalization

{pstd}
{it:Pre-treatment severity} (Section 3.1) quantifies the magnitude 
of parallel trend violations using the normalized Lp-norm:

{p 8 8 2}
S_pre = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p )^{1/p}

{pstd}
{it:Pre-test} (Section 4.2, Theorem 1): Following standard hypothesis testing
convention where phi = 1 indicates rejection of the null hypothesis H0: S_pre <= M,
the test statistic is defined as:

{p 8 8 2}
phi = 1{S_hat_pre > M}    (phi = 0: PASS, phi = 1: FAIL)

{pstd}
The stored result {cmd:e(pretest_pass)} = 1 - phi provides a more intuitive indicator
where 1 indicates that the pre-test passed and extrapolation is justified.

{pstd}
{it:Conditional Extrapolation Assumption} (Assumption 3, Section 3.1) formalizes 
when extrapolation from pre- to post-treatment is justified:

{p 8 8 2}
If S_pre <= M, then S_post <= S_pre

{pstd}
The threshold M reflects the analyst's judgment about the maximum severity of 
pre-treatment violations under which extrapolation remains substantively plausible.

{pstd}
Under this assumption, when the pre-test passes (phi=0), Proposition 1 
establishes a sharp upper bound on the bias: |tau_bar - delta_bar| <= kappa * S_pre 
(iterative mode) or |tau_bar - delta_bar| <= S_pre (overall mode). The constant 
kappa depends on T_post and p through the Holder conjugate relationship.

{pstd}
{bf:Two Modes:}

{pstd}
{it:Iterative mode} (default): Uses period-to-period violations nu_t. The 
confidence interval includes a kappa multiplier that accounts for how violations
can accumulate over post-treatment periods.

{pstd}
{it:Overall mode} ({opt overall} option): Uses cumulative violations nu_bar_t, 
corresponding to TWFE lead coefficients. The confidence interval has no kappa 
multiplier because cumulative violations directly bound the bias.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt treatment(varname)} specifies the binary treatment indicator variable.
Must contain only values 0 (control) and 1 (treated). For panel data, this
can be either time-invariant (group indicator D_i) or time-varying (D_it = 1
when i is treated at time t). Violations of the 0/1 restriction cause 
error 103.

{phang}
{opt time(varname)} specifies the time variable. Must be numeric. If gaps 
exist in the time values, the command creates an internal mapping to 
consecutive indices.

{phang}
{opt threshold(#)} specifies the acceptable violation threshold M > 0. This 
parameter is central to the conditional extrapolation assumption (Assumption 3): 
if S_pre <= M, then post-treatment violations are assumed bounded by pre-treatment 
violations. The choice of M should reflect domain knowledge about the maximum 
magnitude of pre-treatment parallel trend violations under which extrapolation 
remains substantively plausible. Common approaches include: (i) a fraction of 
the outcome variable's standard deviation, (ii) a fraction of the expected 
treatment effect magnitude, or (iii) economic reasoning about acceptable bias. 
Violations of M > 0 cause error 105.

{dlgtab:Optional}

{phang}
{opt treat_time(#)} specifies the treatment time point t0, defined as the
first post-treatment period. Auto-detection is only possible when the treatment
variable is time-varying (i.e., switches from 0 to 1 at time t0). If treatment
is a time-invariant group indicator (the common case in DID), t0 cannot be
inferred from data and this option must be specified explicitly.
Requires T_pre >= 2 (at least 2 pre-treatment periods); violations cause 
error 104.

{phang}
{opt p(#)} specifies the severity norm index p >= 1. Default is p=2 (Euclidean 
norm, standard in econometrics). The choice of p reflects how violations are 
aggregated across pre-treatment periods (Section 3.1): p=1 averages absolute 
violations (robust to outliers); p=2 gives root-mean-square violations (sensitive 
to large deviations); large p approximates the maximum absolute violation 
(worst-case single-period measure). Violations of p >= 1 cause error 106.

{pmore}
{bf:Approximating p=∞ (L-infinity norm):} The paper discusses p=∞ which captures
the maximum absolute violation. Since Stata does not support infinity as a numeric
value, use {cmd:p(1e10)} or {cmd:p(1e6)} as a practical approximation. For example:
{cmd:pretest y, ... p(1e10)} will compute severity as approximately the maximum
absolute violation across pre-treatment periods.

{phang}
{opt alpha(#)} specifies the significance level for confidence intervals.
Must be in the interval (0, 1). Default is 0.05 for 95% confidence intervals.
Takes precedence over {opt level()} if both are specified. Violations cause 
error 107.

{phang}
{opt level(#)} specifies the confidence level as a percentage. Default is 95.
Equivalent to alpha = 1 - level/100. If {opt alpha()} is also specified, 
{opt alpha()} takes precedence.

{phang}
{opt cluster(varname)} specifies a variable identifying clusters for 
cluster-robust standard errors. Typically the panel identifier (e.g., state,
firm) for panel data.

{pmore}
{bf:Understanding Cluster SE Behavior in DID:} Unlike cross-sectional 
regression, the relationship between cluster and non-cluster standard errors 
in DID can be counter-intuitive. Because DID involves time-differencing, 
{bf:time-invariant cluster effects are automatically eliminated}. This has 
important implications:

{pmore}
{bf:When cluster SE < non-cluster SE:} This occurs when cluster effects are 
primarily time-invariant. The non-cluster estimator incorrectly includes the 
variance of these eliminated effects, while the cluster estimator correctly 
reflects only the residual noise variance. This is {bf:expected behavior} consistent 
with the econometric theory of differencing estimators.

{pmore}
{bf:When cluster SE > non-cluster SE:} This occurs when there are time-varying 
cluster effects (cluster-by-period shocks) or within-cluster serial correlation 
that is not fully eliminated by differencing.

{pmore}
{bf:Comparison with Stata regress:} If you compare {cmd:pretest} cluster SE 
with {cmd:regress ... , vce(cluster)} on cross-sectional data (e.g., only 
post-treatment periods), the ratios may differ because {cmd:regress} does not 
involve time-differencing and thus does not eliminate time-invariant effects.
For a proper comparison, use {cmd:reghdfe ... , absorb(id) vce(cluster var)} 
which includes unit fixed effects.

{phang}
{opt overall} specifies that overall (cumulative) violations should be used
instead of iterative (period-to-period) violations. In overall mode, the 
severity measure uses nu_bar_t instead of nu_t, and the confidence interval
has no kappa multiplier. This mode corresponds to testing lead coefficients 
in TWFE event study specifications.

{phang}
{opt nograph} suppresses the event study graph that is displayed by default.

{phang}
{opt simulate(#)} specifies the number of Monte Carlo simulations used to
compute the critical value f(alpha, Sigma). Default is 5000. Larger values
give more precise critical values but take longer to compute. Must be >= 100;
violations cause error 110.

{phang}
{opt seed(#)} specifies the random number seed for Monte Carlo simulations,
ensuring reproducibility. Default is 12345.

{phang}
{opt diagnose} displays detailed diagnostic information after the main output,
including:

{phang2}- Parameter vector theta_hat with all T-1 components{p_end}
{phang2}- Pre-treatment violations nu_hat_t for each period{p_end}
{phang2}- Post-treatment DID estimates delta_hat_t for each period{p_end}
{phang2}- Eigenvalues of the covariance matrix Sigma_hat{p_end}
{phang2}- Monte Carlo simulation parameters and critical value{p_end}

{pmore}
This option is useful for understanding the components that contribute to the
severity measure and for diagnosing potential numerical issues.

{dlgtab:Graph customization}

{phang}
{opt ci_opt_pass(string)} specifies additional or override options for the 
confidence interval lines when the pre-test passes. These options are appended 
to the default styling ({cmd:lcolor(navy)} for pre-treatment, {cmd:lcolor(maroon)} 
for post-treatment), allowing you to override colors, line patterns, and widths.
Example: {cmd:ci_opt_pass(lcolor(black) lwidth(thick))}

{phang}
{opt ci_opt_fail(string)} specifies options for confidence interval lines when 
the pre-test fails. Default styling uses dashed lines with reduced opacity to 
indicate that inference may be invalid. Example: {cmd:ci_opt_fail(lcolor(gs8) lpattern(dot))}

{phang}
{opt line_opt_m(string)} specifies options for the horizontal threshold M lines.
Default styling is {cmd:lcolor(orange%70) lpattern(shortdash)}. These lines mark
the acceptable violation bounds +M and -M on the y-axis.
Example: {cmd:line_opt_m(lcolor(red) lwidth(thick))}

{phang}
{opt marker_opt_pre(string)} specifies options for pre-treatment point estimate 
markers. Default styling is {cmd:mcolor(navy) msymbol(O) msize(medlarge)}.
Example: {cmd:marker_opt_pre(mcolor(black) msymbol(D))}

{phang}
{opt marker_opt_post(string)} specifies options for post-treatment point estimate
markers. Default styling depends on pretest result: solid markers when passing,
hollow markers when failing. Example: {cmd:marker_opt_post(mcolor(green) msymbol(T))}

{phang}
{it:twoway_options} are any of the options documented in {help twoway_options}.
These are passed through to the final {cmd:twoway} graph command and allow full
customization of titles, axes, legends, schemes, and other graph elements.
Common options include:

{phang2}- {cmd:title()}, {cmd:xtitle()}, {cmd:ytitle()} - Customize titles{p_end}
{phang2}- {cmd:xlabel()}, {cmd:ylabel()} - Customize axis labels{p_end}
{phang2}- {cmd:legend()} - Customize legend appearance and position{p_end}
{phang2}- {cmd:scheme()} - Apply a different graph scheme (e.g., {cmd:scheme(s2mono)}){p_end}
{phang2}- {cmd:note()} - Add or suppress footnotes{p_end}

{pmore}
{bf:Override precedence:} User options are appended after default options, so the
Stata rule of "last option wins" ensures user preferences take precedence while
preserving defaults for unspecified options.

{marker examples}{...}
{title:Examples}

{pstd}{bf:Loading Example Data:}{p_end}

{pstd}The package includes the Proposition 99 (California Tobacco Control Program) dataset.
To load it directly from GitHub:{p_end}

{phang2}{cmd:. webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"}{p_end}
{phang2}{cmd:. webuse prop99_smoking.dta, clear}{p_end}

{pstd}{bf:Example 1: Basic panel data analysis}{p_end}

{pstd}Generate simulated panel data and run the pre-test:{p_end}

{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. gen id = ceil(_n/10)}{p_end}
{phang2}{cmd:. gen time = mod(_n-1, 10) + 1}{p_end}
{phang2}{cmd:. gen treat = (id <= 25)}{p_end}
{phang2}{cmd:. gen outcome = rnormal() + treat*(time >= 6)*0.5}{p_end}
{phang2}{cmd:. xtset id time}{p_end}
{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6)}{p_end}

{pstd}Examine stored results:{p_end}

{phang2}{cmd:. ereturn list}{p_end}
{phang2}{cmd:. display "Pre-test passed: " e(pretest_pass)}{p_end}
{phang2}{cmd:. display "Estimated severity: " e(S_pre)}{p_end}

{pstd}{bf:Example 2: Overall violations mode (Less conservative)}{p_end}

{pstd}If the default iterative mode is too strict (due to the kappa multiplier), use the overall violations mode. This mode assumes that cumulative violations are bounded directly, limiting the worst-case bias (kappa = 1) and yielding narrower confidence intervals:{p_end}

{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6) overall}{p_end}

{pstd}Compare with iterative mode - note that overall mode has no kappa multiplier:{p_end}

{phang2}{cmd:. display "Mode: " e(mode)}{p_end}
{phang2}{cmd:. display "Kappa: " e(kappa)}{p_end}

{pstd}{bf:Example 3: Cluster-robust standard errors}{p_end}

{pstd}With clustering at the individual level:{p_end}

{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6) cluster(id)}{p_end}

{pstd}{bf:Example 4: Custom significance level and norm}{p_end}

{pstd}Use 90% confidence level with L1 norm:{p_end}

{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6) alpha(0.10) p(1)}{p_end}

{pstd}{bf:Example 5: Suppress graph and adjust simulations}{p_end}

{pstd}For faster computation or batch processing:{p_end}

{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6) nograph simulate(1000) seed(42)}{p_end}

{pstd}{bf:Example 6: Custom graph styling for publication}{p_end}

{pstd}Create a black and white graph suitable for print journals:{p_end}

{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6) ///}{p_end}
{phang2}{cmd:    ci_opt_pass(lcolor(black) lwidth(medium)) ///}{p_end}
{phang2}{cmd:    marker_opt_pre(mcolor(black) msymbol(O)) ///}{p_end}
{phang2}{cmd:    marker_opt_post(mcolor(black) msymbol(S)) ///}{p_end}
{phang2}{cmd:    line_opt_m(lcolor(gs5) lpattern(dash)) ///}{p_end}
{phang2}{cmd:    scheme(s2mono) title("Event Study")}{p_end}

{pstd}Suppress titles for LaTeX integration:{p_end}

{phang2}{cmd:. pretest outcome, treatment(treat) time(time) threshold(0.5) treat_time(6) ///}{p_end}
{phang2}{cmd:    title("") note("") legend(off)}{p_end}

{pstd}Export graph to publication formats:{p_end}

{phang2}{cmd:. graph export "figure1.png", replace width(2400)}{p_end}
{phang2}{cmd:. graph export "figure1.pdf", replace}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pretest} stores the following in {cmd:e()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 27 2: Scalars}{p_end}
{synopt:{cmd:e(T)}}total number of time periods{p_end}
{synopt:{cmd:e(t0)}}treatment time point (first post-treatment period){p_end}
{synopt:{cmd:e(T_pre)}}number of pre-treatment periods (= t0 - 1){p_end}
{synopt:{cmd:e(T_post)}}number of post-treatment periods (= T - T_pre){p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n)}}number of observations (alias){p_end}
{synopt:{cmd:e(is_panel)}}1 if panel data, 0 if repeated cross-sections{p_end}
{synopt:{cmd:e(p)}}severity norm used{p_end}
{synopt:{cmd:e(alpha)}}significance level{p_end}
{synopt:{cmd:e(level)}}confidence level (%){p_end}
{synopt:{cmd:e(M)}}threshold M (primary name){p_end}
{synopt:{cmd:e(threshold)}}threshold M (alias){p_end}
{synopt:{cmd:e(S_pre)}}estimated pre-treatment severity{p_end}
{synopt:{cmd:e(S_pre_se)}}standard error of S_pre (Delta method){p_end}
{synopt:{cmd:e(kappa)}}kappa constant (iterative mode only){p_end}
{synopt:{cmd:e(f_alpha)}}Monte Carlo critical value f(alpha, Sigma){p_end}
{synopt:{cmd:e(phi)}}pre-test result: 0=PASS (S_pre <= M), 1=FAIL (S_pre > M), .=data invalid{p_end}
{synopt:{cmd:e(data_valid)}}data validity indicator: 1=valid, 0=missing periods detected{p_end}
{synopt:{cmd:e(pretest_pass)}}pre-test passed indicator: 1=PASS, 0=FAIL{p_end}
{synopt:{cmd:e(delta_bar)}}average DID change δ̄ relative to t0 (NOT the ATT level, see note){p_end}
{synopt:{cmd:e(ATT)}}alias for delta_bar ({bf:caution}: this is δ̄, not τ){p_end}
{synopt:{cmd:e(se_delta_bar)}}standard error of δ̄{p_end}
{synopt:{cmd:e(ci_lower)}}conditional CI lower bound (if pretest_pass=1){p_end}
{synopt:{cmd:e(ci_upper)}}conditional CI upper bound (if pretest_pass=1){p_end}
{synopt:{cmd:e(ci_conv_lower)}}conventional CI lower bound (assuming parallel trends){p_end}
{synopt:{cmd:e(ci_conv_upper)}}conventional CI upper bound{p_end}
{synopt:{cmd:e(sims)}}number of Monte Carlo simulations{p_end}
{synopt:{cmd:e(seed)}}random number seed used{p_end}

{synoptset 23 tabbed}{...}
{p2col 5 23 27 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}"pretest"{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}outcome variable name (Stata convention){p_end}
{synopt:{cmd:e(outcome)}}outcome variable name (alias){p_end}
{synopt:{cmd:e(treatment)}}treatment variable name{p_end}
{synopt:{cmd:e(time)}}time variable name{p_end}
{synopt:{cmd:e(mode)}}"iterative" or "overall"{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable name (Stata convention){p_end}
{synopt:{cmd:e(cluster)}}cluster variable name (alias){p_end}
{synopt:{cmd:e(title)}}"Conditional Extrapolation Pre-Test"{p_end}

{synoptset 23 tabbed}{...}
{p2col 5 23 27 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector [ATT] (1 x 1){p_end}
{synopt:{cmd:e(V)}}variance matrix [Var(ATT)] (1 x 1){p_end}
{synopt:{cmd:e(nu)}}iterative violations (T_pre-1 x 1){p_end}
{synopt:{cmd:e(delta)}}DID estimates by period (T_post x 1){p_end}
{synopt:{cmd:e(theta)}}full parameter vector (T-1 x 1){p_end}
{synopt:{cmd:e(Sigma)}}asymptotic covariance matrix (T-1 x T-1){p_end}

{pstd}
{bf:Note on Conditional CI:} Conditional confidence intervals {cmd:e(ci_lower)} and {cmd:e(ci_upper)}
are only stored when the pre-test passes ({cmd:e(pretest_pass)} = 1). When the
pre-test fails, these results are not meaningful and are not stored.

{pstd}
{bf:Note on data validity:} When {cmd:e(phi)} is missing (.), the covariance matrix
could not be computed due to missing observations in some time periods.
Common causes include:

{p 8 12 2}(a) Empty cells in some time-treatment combinations{p_end}
{p 8 12 2}(b) Gaps in the time variable (e.g., years 68, 69, 70, 73, 75 with 71, 72, 74 missing){p_end}
{p 8 12 2}(c) Missing values creating incomplete time-group cells{p_end}

{pstd}
When data is invalid: {cmd:e(S_pre)}, {cmd:e(f_alpha)}, {cmd:e(ci_lower)}, {cmd:e(ci_upper)}
will also be missing, and {cmd:e(pretest_pass)} = 0. Use {cmd:e(data_valid)} to check.

{pstd}
{bf:Important: Understanding delta_t vs tau_t}

{pstd}
{bf:Note:} The estimate delta_bar is the DID estimand, not the traditional ATT.

{pstd}
This command estimates delta_t (the DID estimand), NOT tau_t (the ATT). These 
are distinct quantities with an important mathematical relationship.

{pstd}
{bf:Definitions:}

{phang2}{bf:tau_t (ATT at time t)} = E[Y^(1)_t - Y^(0)_t | D=1]{p_end}
{phang2}The traditional Average Treatment Effect on the Treated at time t.{p_end}

{phang2}{bf:delta_t (DID estimand)} = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]{p_end}
{phang2}The change in outcomes relative to the treatment onset t0.{p_end}

{pstd}
{bf:Mathematical Relationship} (Section 2.1, Equation 176-179):

{p 8 8 2}
tau_t = delta_t + nu_bar_t

{pstd}
where nu_bar_t is the overall violation of parallel trends at time t. This 
decomposition shows that the DID estimand equals the ATT only when parallel 
trends hold exactly (nu_bar_t = 0).

{pstd}
{bf:Key Implication - What delta_t Actually Measures:}

{pstd}
Under the assumption that parallel trends hold (nu_bar_t = 0), we have:

{p 8 8 2}
delta_t = tau_t - tau_{t0}

{pstd}
This fundamental identity means delta_t measures the {bf:change in treatment 
effect relative to t0}, not the treatment effect level. This is the intended 
design of the methodology in Mikhaeil and Harshaw (2025).

{pstd}
{bf:Practical Implications:}

{phang2}- {bf:Constant treatment effect}: If tau_t = tau for all t >= t0
(e.g., treatment adds a constant +2.0), then delta_t = tau - tau = 0.
Thus delta_bar is approximately 0 even when the true ATT is substantial.
This is the intended behavior per Mikhaeil and Harshaw (2025).{p_end}

{phang2}- {bf:Time-varying treatment effect}: If tau_{t0} = 2 and tau_{t0+1} = 4, 
then delta_{t0+1} = 4 - 2 = 2 (captures the growth in effect).{p_end}

{phang2}- {bf:delta_{t0} = 0 by construction}: By definition, delta_{t0} = 
(Y_{t0} - Y_{t0}) - (Y_{t0} - Y_{t0}) = 0 identically.{p_end}

{pstd}
{bf:Confidence Interval Interpretation:}

{pstd}
The conditional CI is centered at delta_bar and adjusted by kappa * S_pre to
account for potential parallel trends violations. Under the conditional extrapolation
assumption (Assumption 3), this CI provides valid asymptotic coverage for tau_bar
(the average ATT), not for delta_bar itself.

{pstd}
{bf:When to Use This Method:}

{phang2}- When parallel trends may not hold exactly and robust inference is needed{p_end}
{phang2}- The method bounds worst-case bias and produces conservative CIs{p_end}
{phang2}- When S_pre is small (order n^{-1/2}), CI width approaches the usual rate{p_end}

{pstd}
{bf:For Traditional ATT Estimation:}

{pstd}
For the traditional ATT (level effect), use standard DID estimators:

{phang2}{cmd:. xtdidregress (y) (d), group(id) time(t)}{p_end}
{phang2}{cmd:. regress y i.treat##i.post}{p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{bf:DID Estimator:}

{p 8 8 2}
delta_t = (Y_bar_{t,D=1} - Y_bar_{t0,D=1}) - (Y_bar_{t,D=0} - Y_bar_{t0,D=0})

{p 8 8 2}
delta_bar = (1/T_post) * sum_{t=t0}^{T} delta_t

{pstd}
{bf:Severity Measure:}

{pstd}
For iterative mode:

{p 8 8 2}
S_pre = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p )^{1/p}

{pstd}
For overall mode:

{p 8 8 2}
S_pre^Delta = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_bar_t|^p )^{1/p}

{pstd}
{bf:Kappa Constant} (Section 3.2, Proposition 1):

{pstd}
The bias bound constant kappa captures how iterative violations accumulate over 
post-treatment periods. It is defined as the normalized Lq-norm of time weights:

{p 8 8 2}
kappa = ( (1/T_post) * sum_{t=1}^{T_post} t^q )^{1/q}

{pstd}
where q is the Holder conjugate of p, defined by 1/p + 1/q = 1. Equivalently, 
q = p/(p-1) for p > 1, q = infinity for p = 1, and q = 1 for p = infinity.

{pstd}
Closed-form solutions for common choices of p:

{phang2}p = 1:     kappa = T_post (worst-case accumulation){p_end}
{phang2}p = 2:     kappa = sqrt((T_post+1)(2*T_post+1)/(6*T_post)) (L2 self-dual){p_end}
{phang2}p = inf:   kappa = (T_post+1)/2 (arithmetic mean of time weights){p_end}

{pstd}
{bf:Conditional Confidence Interval - Iterative Mode} (Section 5.1, Theorem 2):

{pstd}
The confidence interval accounts for both sampling uncertainty and worst-case 
bias under the conditional extrapolation assumption:

{p 8 8 2}
I = delta_bar_hat +/- {kappa * S_pre_hat + f(alpha, Sigma_hat) / sqrt(n)}

{pstd}
where f(alpha, Sigma) is the (1-alpha) quantile of the distribution of psi(Z) 
for Z ~ N(0, Sigma), computed via Monte Carlo simulation.

{pstd}
{bf:Conditional Confidence Interval - Overall Mode} (Appendix C):

{p 8 8 2}
I^Delta = delta_bar_hat +/- {S_pre_hat^Delta + f^Delta(alpha, Sigma_hat^Delta) / sqrt(n)}

{pstd}
{bf:Key difference:} Overall mode has no kappa multiplier because overall 
violations nu_bar_t directly measure cumulative deviation from parallel trends, 
so no additional accumulation factor is needed.

{pstd}
{bf:Conditional Coverage Guarantee} (Section 5.1, Theorem 2):

{pstd}
Under the conditional extrapolation assumption (Assumption 3) and well-separated 
null (|S_pre - M| = omega(n^{-1/2})), the confidence interval achieves:

{p 8 8 2}
liminf_{n->infinity} Pr{tau_bar in I | phi=0} >= 1 - alpha

{pstd}
This conditional validity property resolves the coverage distortions identified 
by Roth (2022), where conventional CIs suffer from severe undercoverage when 
conditioned on passing preliminary tests for parallel trends.


{marker references}{...}
{title:References}

{phang}
Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
Inference when Testing Violations of Parallel Trends for Difference-in-
Differences. {it:arXiv preprint} arXiv:2510.26470.
{browse "https://arxiv.org/abs/2510.26470"}

{phang}
Rambachan, A. and J. Roth. 2023. A More Credible Approach to Parallel
Trends. {it:Review of Economic Studies} 90(5): 2555-2591.
{browse "https://doi.org/10.1093/restud/rdad018"}

{phang}
Roth, J. 2022. Pre-test with Caution: Event-Study Estimates after Testing
for Parallel Trends. {it:American Economic Review: Insights} 4(3): 305-322.
{browse "https://doi.org/10.1257/aeri.20210236"}


{marker author}{...}
{title:Authors}

{pstd}
{bf:Stata Implementation:}

{phang2}Xuanyu Cai, City University of Macau{p_end}
{phang2}Email: {browse "mailto:xuanyuCAI@outlook.com":xuanyuCAI@outlook.com}{p_end}

{phang2}Wenli Xu, City University of Macau{p_end}
{phang2}Email: {browse "mailto:wlxu@cityu.edu.mo":wlxu@cityu.edu.mo}{p_end}

{pstd}
{bf:Methodology:} Jonas M. Mikhaeil and Christopher Harshaw,
Department of Statistics, Columbia University.

{pstd}
Please cite both the methodology paper and this Stata implementation when using
this command in published work.
