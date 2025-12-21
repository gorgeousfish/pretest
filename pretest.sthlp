{smcl}
{* *! version 0.1.0  21dec2025}{...}
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
{bf:pretest} {hline 2} Conditional extrapolation pre-test for difference-in-differences designs


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:pretest}
{it:outcome}
{cmd:,}
{opt treat:ment(varname)}
{opt time(varname)}
{opt thr:eshold(#)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt treat:ment(varname)}}binary treatment indicator (0/1){p_end}
{synopt:{opt time(varname)}}time variable{p_end}
{synopt:{opt thr:eshold(#)}}acceptable violation threshold M > 0{p_end}

{syntab:Optional}
{synopt:{opt treat_time(#)}}treatment time point t0; auto-detected only if treatment is time-varying{p_end}
{synopt:{opt p(#)}}severity norm p >= 1; default is {cmd:p(2)}{p_end}
{synopt:{opt alpha(#)}}significance level; default is {cmd:alpha(0.05)}{p_end}
{synopt:{opt level(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt cl:uster(varname)}}cluster variable for standard errors{p_end}
{synopt:{opt over:all}}use overall violations mode (no kappa multiplier){p_end}
{synopt:{opt nog:raph}}suppress event study graph{p_end}
{synopt:{opt sim:ulate(#)}}Monte Carlo simulations; default is {cmd:simulate(5000)}{p_end}
{synopt:{opt seed(#)}}random number seed; default is {cmd:seed(12345)}{p_end}
{synopt:{opt diag:nose}}display detailed diagnostic information{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{bf:{err:IMPORTANT REQUIREMENTS:}}

{phang2}{err:•} {bf:Minimum 3 time periods required}: T_pre >= 2 (at least 2 pre-treatment periods){p_end}
{phang2}{err:•} {bf:Block adoption design}: All treated units must receive treatment at the same time t0{p_end}
{phang2}{err:•} {bf:Binary treatment}: Treatment variable must be 0/1{p_end}

{pstd}
{bf:Note:} For datasets with only 2 time periods (e.g., pre vs post), this command
cannot be used because iterative violations ν_t are only defined for t >= 2.
In such cases, use standard DID estimators instead:

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
{it:Iterative violation} at time t measures the period-to-period deviation from 
parallel trends:

{p 8 8 2}
nu_t = E[Y(0)_t - Y(0)_{t-1} | D=1] - E[Y(0)_t - Y(0)_{t-1} | D=0]

{pstd}
{it:Overall violation} at time t is the cumulative sum of iterative violations,
equivalent to TWFE event study lead coefficients:

{p 8 8 2}
nu_bar_t = sum_{s=2}^{t} nu_s    (with nu_bar_1 = 0)

{pstd}
{it:Pre-treatment severity} measures the magnitude of violations using the Lp norm:

{p 8 8 2}
S_pre = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p )^{1/p}

{pstd}
{it:Pre-test} (Theorem 1). {bf:Note:} This implementation uses
phi = 1{S_pre > M}, where phi=1 means FAIL (reject extrapolation) and phi=0 means PASS:

{p 8 8 2}
phi = 1{S_pre > M}    (phi=0: PASS, phi=1: FAIL)

{pstd}
The returned {cmd:e(pretest_pass)} = 1 - phi provides a more intuitive indicator
where 1 means the pre-test passed.

{pstd}
{it:Conditional Extrapolation Assumption} (Assumption 3):

{p 8 8 2}
If S_pre <= M, then S_post <= S_pre

{pstd}
Under this assumption, if the pre-test passes, the worst-case bias in the DID 
estimator is bounded by kappa * S_pre (for iterative mode) or S_pre (for overall 
mode), where kappa depends on T_post and p.

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
{opt threshold(#)} specifies the acceptable violation threshold M. Must be 
positive (M > 0). If S_pre > M, the pre-test fails and confidence intervals
are not reported. The threshold should be chosen based on domain knowledge 
about what magnitude of parallel trend violations would be acceptable for
your application. A common approach is to use a fraction of the outcome
variable's scale. Violations of M > 0 cause error 105.

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
{opt p(#)} specifies the severity norm p >= 1. Default is p=2 (Euclidean/L2 norm).
Use p=1 for average absolute violations (L1 norm). Use large values (e.g., p=100)
to approximate the maximum absolute violation (L-infinity norm). Violations of 
p >= 1 cause error 106.

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
reflects only the residual noise variance. This is {bf:correct behavior}, not 
a bug.

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

{phang2}• Parameter vector θ̂ with all T-1 components{p_end}
{phang2}• Pre-treatment violations ν̂_t for each period{p_end}
{phang2}• Post-treatment DID estimates δ̂_t for each period{p_end}
{phang2}• Eigenvalues of the covariance matrix Σ̂{p_end}
{phang2}• Monte Carlo simulation parameters and critical value{p_end}

{pmore}
This option is useful for understanding the components that contribute to the
severity measure and for diagnosing potential numerical issues.


{marker examples}{...}
{title:Examples}

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

{pstd}{bf:Example 2: Overall violations mode}{p_end}

{pstd}Use overall violations mode (equivalent to testing TWFE lead coefficients):{p_end}

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
{synopt:{cmd:e(phi)}}pre-test result: 0=PASS, 1=FAIL, 2=DATA_INVALID{p_end}
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
{bf:{err:Note on phi=2 (DATA_INVALID):}} When {cmd:e(phi)} = 2, this indicates that the
covariance matrix could not be computed due to missing observations in some time periods.
This typically occurs when:

{p 8 12 2}(a) Some time periods have no observations for the treatment or control group{p_end}
{p 8 12 2}(b) The data has gaps in the time variable (e.g., years 68,69,70,73,75 with 71,72,74 missing){p_end}
{p 8 12 2}(c) Missing values in the treatment or outcome variables create empty cells{p_end}

{pstd}
When {cmd:e(phi)} = 2: {cmd:e(S_pre)}, {cmd:e(f_alpha)}, {cmd:e(ci_lower)}, {cmd:e(ci_upper)} will be missing (.),
and {cmd:e(pretest_pass)} = 0. Use {cmd:e(data_valid)} to explicitly check data validity.

{pstd}
{bf:{err:CRITICAL: Understanding δ_t vs τ_t - READ BEFORE INTERPRETING RESULTS}}

{pstd}
{err:WARNING: The estimate δ̄ (delta_bar) is NOT the traditional ATT!}

{pstd}
This command estimates δ_t (the DID estimand), NOT τ_t (the ATT). These are 
different quantities with an important mathematical relationship.

{pstd}
{bf:Definitions:}

{phang2}{bf:τ_t (ATT at time t)} = E[Y^(1)_t - Y^(0)_t | D=1]{p_end}
{phang2}This is the traditional Average Treatment Effect on the Treated at time t.{p_end}

{phang2}{bf:δ_t (DID estimand)} = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]{p_end}
{phang2}This measures the change in outcomes relative to t0.{p_end}

{pstd}
{bf:Mathematical Relationship:}

{p 8 8 2}
τ_t = δ_t + ν̄_t

{pstd}
where ν̄_t is the overall violation of parallel trends at time t.

{pstd}
{bf:Key Implication - What δ_t Actually Measures:}

{pstd}
Under the assumption that parallel trends hold (ν̄_t = 0), we have:

{p 8 8 2}
δ_t = τ_t - τ_{t0}

{pstd}
This means δ_t measures the {bf:change in treatment effect relative to t0}, not 
the treatment effect itself!

{pstd}
{bf:Practical Consequences:}

{phang2}• {bf:Constant treatment effect (IMPORTANT!)}: If τ_t = τ for all t ≥ t0 
(e.g., treatment always adds +2.0), then δ_t = τ - τ = 0. 
{bf:This means δ̄ ≈ 0 even when the true ATT is large!} This is correct behavior 
based on the paper's definition, not a bug.{p_end}

{phang2}• {bf:Growing treatment effect}: If τ_4 = 2 and τ_5 = 4, then δ_5 = 4 - 2 = 2.{p_end}

{phang2}• {bf:δ_{t0} = 0 always}: By definition, δ_{t0} = (Y_{t0} - Y_{t0}) - (Y_{t0} - Y_{t0}) = 0.{p_end}

{pstd}
{bf:Confidence Interval Interpretation:}

{pstd}
The conditional CI from this command is centered at δ̄ and adjusted by κ·S_pre to 
account for potential parallel trends violations. Under the conditional extrapolation
assumption, this CI has valid coverage for τ̄ (the average ATT), NOT for δ̄.

{pstd}
{bf:When to Use This Method:}

{phang2}• Best suited when you suspect parallel trends may not hold exactly{p_end}
{phang2}• The method bounds the bias from violations and produces conservative CIs{p_end}
{phang2}• When S_pre ≈ 0 (perfect parallel trends), CI width is minimized{p_end}

{pstd}
{bf:For Traditional ATT Estimation:}

{pstd}
If you want the traditional ATT (level effect), use standard DID estimators:

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
{bf:Kappa Constant} (Proposition 1):

{p 8 8 2}
kappa = ( (1/T_post) * sum_{t=1}^{T_post} t^q )^{1/q}, where 1/p + 1/q = 1

{pstd}
Special cases: p=1 gives kappa=T_post; p=2 gives kappa=sqrt((T_post+1)(2T_post+1)/6);
p=infinity gives kappa=(T_post+1)/2.

{pstd}
{bf:Conditional Confidence Interval - Iterative Mode} (Theorem 2):

{p 8 8 2}
CI = delta_bar +/- {kappa * S_pre + f(alpha, Sigma) / sqrt(n)}

{pstd}
{bf:Conditional Confidence Interval - Overall Mode:}

{p 8 8 2}
CI^Delta = delta_bar +/- {S_pre^Delta + f^Delta(alpha, Sigma^Delta) / sqrt(n)}

{pstd}
{bf:Key difference:} Overall mode has no kappa multiplier because cumulative 
violations directly enter the bias bound.

{pstd}
{bf:Conditional Coverage Guarantee} (Theorem 2):

{p 8 8 2}
liminf_{n->infinity} Pr{tau_bar in CI | phi=0} >= 1 - alpha

{pstd}
This resolves the conditional coverage problem identified by Roth (2022).


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

{phang}
Roth, J. 2022. Pre-test with Caution: Event-Study Estimates after Testing
for Parallel Trends. {it:American Economic Review: Insights} 4(3): 305-322.


{marker author}{...}
{title:Authors}

{pstd}
{bf:Stata implementation:} Xuanyu Cai and Wenli Xu

{pstd}
{bf:Methodology:} Jonas M. Mikhaeil and Christopher Harshaw,
Department of Statistics, Columbia University.

{pstd}
Please cite both the methodology paper and the Stata implementation when using
this command in published work.


{title:Also see}

{psee}
Manual: {manlink R regress}, {manlink XT xtreg}, {manlink XT xtset}

{psee}
{space 2}Help: {help regress}, {help xtreg}, {help xtset}, {help xtdidregress}
{p_end}
