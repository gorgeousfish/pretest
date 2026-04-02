{smcl}
{* *! version 0.1.1  02apr2026}{...}
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
{viewerjumpto "Citation" "pretest##citation"}{...}

{title:Title}

{phang}
{bf:pretest} {hline 2} Stata module for thresholded DID diagnostics under conditional extrapolation


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:pretest}
{it:outcome}
[{it:if}]
[{it:in}]
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
{synopt:{opt treat_time(#)}}common treatment-onset time t0; auto-detected only from a common 0-to-1 switch path{p_end}
{synopt:{opt p(#)}}severity norm p >= 1; default is {cmd:p(2)}{p_end}
{synopt:{opt al:pha(#)}}nominal significance level for interval construction; if omitted, follows Stata's current confidence level{p_end}
{synopt:{opt level(#)}}nominal confidence level (%) for interval construction; if omitted, follows Stata's current confidence level{p_end}
{synopt:{opt cl:uster(varname)}}cluster variable for robust standard errors after explicit sample prep{p_end}
{synopt:{opt over:all}}use cumulative-violation path (no kappa multiplier){p_end}
{synopt:{opt nog:raph}}suppress event study graph{p_end}
{synopt:{opt sim:ulate(#)}}Monte Carlo simulations; default is {cmd:simulate(5000)}{p_end}
{synopt:{opt seed(#)}}random number seed; default is {cmd:seed(12345)}{p_end}
{synopt:{opt diag:nose}}display detailed diagnostic information{p_end}

{syntab:Graph customization}
{synopt:{opt ci_opt_pass(string)}}override styling for pre-treatment CIs and pass-case post-treatment CIs{p_end}
{synopt:{opt ci_opt_fail(string)}}override CI style in non-pass cases{p_end}
{synopt:{opt line_opt_m(string)}}override threshold M line style{p_end}
{synopt:{opt marker_opt_pre(string)}}override pre-treatment marker style{p_end}
{synopt:{opt marker_opt_post(string)}}override post-treatment marker style{p_end}
{synopt:{it:twoway_options}}any standard Stata {help twoway_options}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:pretest} asks whether the observed pre-treatment path is mild enough to
carry {cmd:delta_bar} into conditional extrapolation in
common-treatment-time block-adoption designs.

{pstd}
One run measures pre-treatment severity, stores a pass indicator, reports the
DID summary {cmd:delta_bar}, and, when {cmd:e(pretest_pass)} = 1, makes
pass-case interval endpoints available to be read as centered on
{cmd:delta_bar} and as informative about the average ATT only under the
conditional extrapolation assumption.

{pstd}
In the shipped Prop99 example, lowering M from 5 to 1 leaves {cmd:delta_bar}
and {cmd:S_pre} unchanged but moves {cmd:e(pretest_pass)} from 1 to 0. The
threshold changes what the user may claim before it changes the DID summary.

{pstd}
In practice, the command separates three choices that applied DID work often
blurs: threshold choice, sample choice, and violation-summary choice.

{pstd}
Read one run in four steps: {cmd:S_pre}, the pass indicator, {cmd:delta_bar},
and any pass-case interval endpoints. When those endpoints appear, read them
as centered on {cmd:delta_bar}. Treat them as informative about the average
ATT only under the conditional extrapolation assumption and the asymptotic
conditions stated in the method paper.

{pstd}
Users choose threshold M under the {it:conditional extrapolation assumption},
so post-treatment interpretation remains conditional on pre-treatment
violations staying below that bound.

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
{it:Pre-test indicator}: At the interface level, use the coding

{p 8 8 2}
phi = 1{S_hat_pre > M}    (phi = 0: pass, phi = 1: fail)

{pstd}
The stored result {cmd:e(pretest_pass)} is the main reporting flag: it equals
1 when phi = 0 and 0 in all non-pass cases. When {cmd:e(pretest_pass)} = 0,
use {cmd:e(phi)} together with {cmd:e(data_valid)} to distinguish
threshold-fail from invalid-data cases.

{pstd}
{it:Conditional Extrapolation Assumption} (Assumption 3, Section 3.1) formalizes 
when extrapolation from pre- to post-treatment is justified:

{p 8 8 2}
If S_pre <= M, then S_post <= S_pre

{pstd}
The threshold M reflects the analyst's judgment about the maximum severity of 
pre-treatment violations under which extrapolation remains substantively plausible.

{pstd}
Under this assumption, when the pre-test passes (phi=0), the method bounds the
bias between tau_bar and delta_bar. In iterative mode, that bound uses
kappa * S_pre. In overall mode, it uses the cumulative-violation object
S_pre^Delta rather than the iterative path with kappa simply turned off. The
constant kappa depends on T_post and p through the Holder conjugate relationship.

{pstd}
{bf:Two Modes:}

{pstd}
{it:Iterative mode} (default): Uses period-to-period violations nu_t. The 
confidence interval includes a kappa multiplier that accounts for how violations
can accumulate over post-treatment periods.

{pstd}
{it:Overall mode} ({opt overall} option): Uses cumulative violations nu_bar_t,
corresponding to TWFE lead coefficients. Its bias bound and interval
construction follow the cumulative-violation path based on S_pre^Delta rather
than the iterative kappa * S_pre path. Use this mode when cumulative drift is
the empirical question.


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
Violations of M > 0 trigger a parser error.

{dlgtab:Optional}

{phang}
{opt treat_time(#)} specifies the treatment time point t0, defined as the
first post-treatment period. Auto-detection is only possible when the
treatment path exhibits one common 0-to-1 switch date. If treatment is a
time-invariant group indicator, or if units do not share one common switch
date, t0 cannot be inferred from data and this option must be specified
explicitly.
Requires T_pre >= 2 (at least 2 pre-treatment periods); violations cause 
error 104.

{phang}
{opt p(#)} specifies the severity norm index p >= 1. Default is p=2 (Euclidean 
norm, standard in econometrics). The choice of p reflects how violations are 
aggregated across pre-treatment periods (Section 3.1): p=1 averages absolute 
violations (robust to outliers); p=2 gives root-mean-square violations (sensitive 
to large deviations); large p approximates the maximum absolute violation 
(worst-case single-period measure). Violations of p >= 1 trigger a parser
error.

{pmore}
{bf:Approximating p=∞ (L-infinity norm):} The paper discusses p=∞ which captures
the maximum absolute violation. Since Stata does not support infinity as a numeric
value, use {cmd:p(1e10)} or {cmd:p(1e6)} as a practical approximation. For example:
{cmd:pretest y, ... p(1e10)} will compute severity as approximately the maximum
absolute violation across pre-treatment periods.

{phang}
{opt alpha(#)} specifies the nominal significance level used in interval
construction. Must be in the interval (0, 1). If {opt alpha()} and
{opt level()} are both omitted, the command inherits Stata's current
confidence level and uses the implied nominal significance level. If both are
specified, {opt alpha()} takes precedence. Violations trigger a parser error.

{phang}
{opt level(#)} specifies the nominal confidence level as a percentage for
interval construction. If {opt alpha()} and {opt level()} are both omitted,
the command inherits Stata's current confidence level. Equivalent to
alpha = 1 - level/100. If {opt alpha()} is also specified,
{opt alpha()} takes precedence.

{phang}
{opt cluster(varname)} specifies a variable identifying clusters for 
cluster-robust standard errors. Typically the panel identifier (e.g., state,
firm) for panel data.

{pmore}
Use clustered conditional inference only after explicitly deleting rows with
missing outcome, treatment, or time values. In the current workflow, leaving
that cleanup to the command's internal complete-case filter can change
covariance-dependent objects such as {cmd:e(Sigma)}, {cmd:f_alpha}, and the
pass-case interval even when {cmd:e(delta_bar)} and {cmd:e(S_pre)} stay
fixed. The safer practice is to prepare the analysis sample first and then
compare clustered and non-clustered runs on that same sample.

{pmore}
Before clustered runs, subset the analysis sample explicitly and rerun
{cmd:pretest} on that fixed sample; inline {cmd:[if] [in]} restrictions are
not a supported substitute for this workflow.

{phang}
{opt overall} specifies the cumulative-violation path instead of iterative
(period-to-period) violations. In overall mode, the severity measure uses
nu_bar_t instead of nu_t, and pass-case interval construction has no kappa
multiplier. This mode corresponds to a cumulative-violation summary related
to TWFE lead coefficients.

{phang}
{opt nograph} suppresses the event study graph that is displayed by default.

{phang}
{opt simulate(#)} specifies the number of Monte Carlo simulations used to
compute the critical value f(alpha, Sigma). Default is 5000. Larger values
give more precise critical values but take longer to compute. Must be >= 100;
violations trigger a parser error.

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
{opt ci_opt_pass(string)} specifies options for the CI layers managed by the
pass-style path. The live graph helper applies this override to the
pre-treatment CI layer and, when the run passes, to the post-treatment
pass-case endpoints. When specified, this completely replaces the default
styling for those layers. You must provide all desired styling options.
Example: {cmd:ci_opt_pass(lcolor(black) lwidth(thick))}

{phang}
{opt ci_opt_fail(string)} specifies options for confidence interval lines in 
non-pass cases. Default styling uses dashed lines with reduced opacity to 
indicate that pass-case conditional inference is unavailable. Example:
{cmd:ci_opt_fail(lcolor(gs8) lpattern(dot))}

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
markers. Default styling depends on run status: solid markers in pass cases and
hollow markers in non-pass cases. Example:
{cmd:marker_opt_post(mcolor(green) msymbol(T))}

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
{bf:Override behavior:} Element-specific options ({cmd:ci_opt_pass}, {cmd:marker_opt_pre}, etc.) 
completely replace the default styling when specified. Unspecified elements retain 
their defaults. Standard {it:twoway_options} are appended and follow Stata's 
"last option wins" rule.

{marker examples}{...}
{title:Examples}

{pstd}{bf:Loading Example Data:}{p_end}

{pstd}The package includes the Proposition 99 (California Tobacco Control Program) dataset.
If {cmd:prop99_smoking.dta} is available locally, either from a repo checkout or
from {cmd:net get pretest}, load that local copy first:{p_end}

{phang2}{cmd:. capture noisily findfile prop99_smoking.dta}{p_end}
{phang2}{cmd:. if _rc == 0 use "`r(fn)'", clear}{p_end}

{pstd}If no local copy is available, fall back to GitHub:{p_end}

{phang2}{cmd:. webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"}{p_end}
{phang2}{cmd:. webuse prop99_smoking.dta, clear}{p_end}

{pstd}A clean {cmd:net install pretest} alone does not place discoverable example
assets on the installed adopath; {cmd:net get pretest} downloads the example
files into the current working directory as loose files.{p_end}

{pstd}{bf:Worked Prop99 illustration}{p_end}

{pstd}The shipped Prop99 script teaches one workflow in four steps. The first
two runs keep the full sample fixed and isolate the reporting gate. The third
run trims the panel and restarts the question on a new sample. The fourth run
keeps that trimmed panel fixed and asks the {cmd:overall} version of the same
question. Clustered reruns stay outside that ladder as fixed-sample
preparation, not as a fifth worked example.{p_end}

{phang2}{cmd:. xtset state year}{p_end}
{phang2}{cmd:. pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)}{p_end}
{phang2}{cmd:. pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(1) nograph}{p_end}
{phang2}{cmd:. keep if year >= 1985 & year <= 1995}{p_end}
{phang2}{cmd:. pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)}{p_end}
{phang2}{cmd:. pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5) overall}{p_end}

{pstd}{bf:Step 1: Baseline run on the full panel.} With {cmd:threshold(5)},
the command returns {cmd:N = 1209}, {cmd:delta_bar = -14.4450},
{cmd:S_pre = 2.1959}, and {cmd:pretest_pass = 1}. Under that gate, users may
discuss the pass-case endpoints, but the center of that discussion remains
{cmd:delta_bar}.{p_end}

{pstd}{bf:Step 2: Tighten the threshold on the same sample.} Lowering the
threshold from 5 to 1 leaves {cmd:delta_bar} and {cmd:S_pre} unchanged but
moves {cmd:pretest_pass} from 1 to 0, so {cmd:e(ci_lower)} and
{cmd:e(ci_upper)} disappear. The threshold changes what the user may claim
before it changes the DID summary.{p_end}

{pstd}{bf:Step 3: Trim the window, then restart the question.} After
{cmd:keep if year >= 1985 & year <= 1995}, the rerun returns
{cmd:N = 429}, {cmd:delta_bar = -9.3421}, and {cmd:S_pre = 2.1767}. That is a
new sample-level question, not the old sample with one option changed. Subset
the data first rather than treat inline {cmd:[if] [in]} as a substitute for a
fixed-sample rerun.{p_end}

{pstd}{bf:Step 4: Ask the overall-mode version on that trimmed sample.} The
final call keeps the 1985--1995 panel fixed and changes only the violation
summary by adding {cmd:overall}. It returns {cmd:N = 429},
{cmd:delta_bar = -9.3421}, {cmd:S_pre = 2.0001}, and a conditional interval of
{cmd:[-25.1525, 6.4683]}. Read this as the cumulative-violation version of the
same trimmed-sample question.{p_end}

{pstd}{bf:Clustered reruns come after the ladder.} First delete rows with
missing outcome, treatment, or time values, then rerun {cmd:pretest} on that
cleaned sample:{p_end}

{phang2}{cmd:. drop if missing(cigsale, treated, year)}{p_end}
{phang2}{cmd:. pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5) cluster(state)}{p_end}

{pstd}That rerun checks clustered inference on a fixed sample. It is not a
second substantive application.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pretest} stores the following in {cmd:e()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 27 2: Scalars}{p_end}
{synopt:{cmd:e(T)}}total number of time periods{p_end}
{synopt:{cmd:e(t0)}}internal treatment-time index (first post-treatment period in consecutive mapped scale; see {cmd:e(t0_orig)} for original calendar value){p_end}
{synopt:{cmd:e(t0_orig)}}treatment time in original calendar units (e.g., the calendar year supplied to {cmd:treat_time()}){p_end}
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
{synopt:{cmd:e(kappa)}}bias-bound scalar kappa; posted as 1 in overall mode{p_end}
{synopt:{cmd:e(f_alpha)}}Monte Carlo critical value f(alpha, Sigma){p_end}
{synopt:{cmd:e(phi)}}pre-test result: 0=PASS (S_pre <= M), 1=FAIL (S_pre > M), .=data invalid{p_end}
{synopt:{cmd:e(data_valid)}}data validity indicator: 1=valid, 0=data issue detected; read with {cmd:e(phi)} and the note below{p_end}
{synopt:{cmd:e(pretest_pass)}}pre-test pass indicator: 1=PASS, 0=NON-PASS; use {cmd:e(phi)} and {cmd:e(data_valid)} to separate fail from invalid data{p_end}
{synopt:{cmd:e(delta_bar)}}average DID change δ̄ relative to t0 (NOT the ATT level, see note){p_end}
{synopt:{cmd:e(ATT)}}compatibility alias mirroring {cmd:e(delta_bar)}; not an ATT-level return{p_end}
{synopt:{cmd:e(se_delta_bar)}}standard error of δ̄{p_end}
{synopt:{cmd:e(ci_lower)}}pass-case interval lower endpoint; read as part of a {cmd:e(delta_bar)}-centered interval when {cmd:e(pretest_pass)} = 1 and as informative about the average ATT only under the conditional extrapolation assumption{p_end}
{synopt:{cmd:e(ci_upper)}}pass-case interval upper endpoint; read as part of a {cmd:e(delta_bar)}-centered interval when {cmd:e(pretest_pass)} = 1 and as informative about the average ATT only under the conditional extrapolation assumption{p_end}
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
{synopt:{cmd:e(b)}}compatibility coefficient vector mirroring {cmd:e(delta_bar)} (1 x 1){p_end}
{synopt:{cmd:e(V)}}compatibility variance matrix for {cmd:e(delta_bar)} (1 x 1){p_end}
{synopt:{cmd:e(nu)}}pre-treatment violation vector (T_pre-1 x 1); iterative nu_t in iterative mode, cumulative nubar_t in overall mode{p_end}
{synopt:{cmd:e(delta)}}DID estimates by period (T_post x 1){p_end}
{synopt:{cmd:e(theta)}}full parameter vector (T-1 x 1); pre-treatment block is nu_t (iterative mode) or nubar_t (overall mode){p_end}
{synopt:{cmd:e(Sigma)}}asymptotic covariance matrix (T-1 x T-1); in overall mode stored in overall-frame coordinates (A*Sigma*A'){p_end}

{pstd}
{bf:Note on Conditional CI:} Conditional confidence intervals {cmd:e(ci_lower)} and {cmd:e(ci_upper)}
are only stored when the pre-test passes ({cmd:e(pretest_pass)} = 1). Read
these pass-case endpoints as centered on {cmd:delta_bar} and as informative
about the average ATT only under the conditional extrapolation assumption. When
the run does not pass, these results are not stored; use {cmd:e(phi)} together
with {cmd:e(data_valid)} to separate threshold-fail from invalid-data cases.

{pstd}
{bf:Note on compatibility returns:} {cmd:e(ATT)}, {cmd:e(b)}, and {cmd:e(V)}
remain available for Stata compatibility. They mirror {cmd:e(delta_bar)} and
its variance, so use {cmd:e(delta_bar)} as the substantive DID summary and do
not read those ATT-labeled surfaces as separate ATT-level returns.

{pstd}
{bf:Note on data validity:} When {cmd:e(phi)} is missing (.), the covariance
matrix could not be computed. Common causes include:

{p 8 12 2}(a) Empty cells in some time-treatment combinations{p_end}
{p 8 12 2}(b) Missing values creating incomplete time-group cells{p_end}
{p 8 12 2}(c) Sample restrictions that remove one group from a period{p_end}

{pstd}
When data is invalid, {cmd:e(S_pre)}, {cmd:e(f_alpha)}, {cmd:e(ci_lower)}, and
{cmd:e(ci_upper)} will also be missing, and {cmd:e(pretest_pass)} = 0. Use
{cmd:e(data_valid)} together with {cmd:e(phi)} to distinguish invalid-data
cases from threshold-fail cases.

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
decomposition isolates the parallel-trend-violation term; delta_t remains the
package's t0-relative DID object.

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
The pass-case interval is centered at delta_bar and adjusted by the relevant
bias term. Under the conditional extrapolation assumption and the asymptotic
conditions stated in the method paper, read this pass-case interval as
informative about tau_bar (the average ATT).

{pstd}
{bf:When to Use This Method:}

{phang2}- When a common-treatment-time block-adoption DID design needs an explicit severity threshold before moving beyond {cmd:delta_bar}{p_end}
{phang2}- When a pass indicator and a pass-case interval reading rule are useful for reporting{p_end}
{phang2}- When the workflow will check severity first, the gate second, and {cmd:delta_bar} before any pass-case interval interpretation{p_end}

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
{bf:Method-paper conditional validity result} (Section 5.1, Theorem 2):

{pstd}
Under the conditional extrapolation assumption (Assumption 3) and well-separated 
null (|S_pre - M| = omega(n^{-1/2})), the method paper establishes:

{p 8 8 2}
liminf_{n->infinity} Pr{tau_bar in I | phi=0} >= 1 - alpha

{pstd}
In the method paper, this conditional validity result addresses the coverage
distortions identified by Roth (2022), where conventional CIs suffer from
severe undercoverage when conditioned on passing preliminary tests for
parallel trends.


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


{marker citation}{...}
{title:Citation}

{pstd}
If you use this package in your research, please cite both the methodology paper
and the Stata implementation:

{pstd}
{bf:APA Format:}

{phang2}Cai, X., & Xu, W. (2025). {it:pretest: Stata module to implement the conditional}
{it:extrapolation pre-test for difference-in-differences} (Version 0.1.1) [Computer software].
GitHub. {browse "https://github.com/gorgeousfish/pretest"}{p_end}

{phang2}Mikhaeil, J. M., & Harshaw, C. (2025). In Defense of the Pre-Test: Valid Inference
when Testing Violations of Parallel Trends for Difference-in-Differences.
{it:arXiv preprint arXiv:2510.26470}. {browse "https://arxiv.org/abs/2510.26470"}{p_end}

{pstd}
{bf:BibTeX:}

        @software{c -(}pretest2025stata,
              title = {c -(}pretest: Stata module to implement the conditional
                       extrapolation pre-test for difference-in-differences{c )-},
              author = {c -(}Xuanyu Cai and Wenli Xu{c )-},
              year = {c -(}2025{c )-},
              version = {c -(}0.1.1{c )-},
              url = {c -(}https://github.com/gorgeousfish/pretest{c )-}
        {c )-}

        @misc{c -(}mikhaeil2025defensepretestvalidinference,
              title = {c -(}In Defense of the Pre-Test: Valid Inference when
                       Testing Violations of Parallel Trends for
                       Difference-in-Differences{c )-},
              author = {c -(}Jonas M. Mikhaeil and Christopher Harshaw{c )-},
              year = {c -(}2025{c )-},
              eprint = {c -(}2510.26470{c )-},
              archivePrefix = {c -(}arXiv{c )-},
              primaryClass = {c -(}stat.ME{c )-},
              url = {c -(}https://arxiv.org/abs/2510.26470{c )-}
        {c )-}
