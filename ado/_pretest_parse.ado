*! _pretest_parse v0.1.0 - Parameter Parsing and Validation Module
*!
*! Description:
*!   Parses and validates pretest command parameters.
*!
*! Syntax:
*!   pretest outcome, treatment(varname) time(varname) threshold(#)
*!
*! Options:
*!   p(#) alpha(#) cluster(varname) nograph simulate(#) overall seed(#)
*!
*! Error codes:
*!   100 - Missing threshold
*!   103 - Treatment not binary 0/1
*!   105 - threshold <= 0
*!   106 - p < 1
*!   107 - alpha not in (0,1)
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

program define _pretest_parse, rclass
    version 17.0

    // Parse command line (supports if/in for validation only)
    // Use -999 as sentinel value to detect if user specified threshold/alpha
    syntax varname(numeric) [if] [in], ///
        TREATment(varname numeric) ///
        TIME(varname numeric) ///
        [THREshold(real -999) ///
         P(real 2) ///
         ALpha(real -999) ///
         CLuster(varname numeric) ///
         TREat_time(real -999) ///
         NOGraph ///
         SIMulate(integer 5000) ///
         OVERall ///
         SEED(integer 12345) ///
         LEVel(cilevel)]

    // Missing threshold: custom error code 100
    if `threshold' == -999 {
        di as error "Error 100: option threshold() required"
        di as error "  Hint: threshold(#) specifies the acceptable level M for pre-treatment violations"
        exit 100
    }

    // alpha/level compatibility handling:
    // If user specified alpha, use alpha
    // If user didn't specify alpha but specified level, compute alpha from level
    // If neither specified, use default 0.05
    if `alpha' == -999 {
        // User did not specify alpha
        if "`level'" != "" {
            // User specified level, compute alpha from level
            local alpha = 1 - `level'/100
        }
        else {
            // Neither specified, use default
            local alpha = 0.05
        }
    }

    // Parameter validation
    if `threshold' <= 0 {
        di as error "threshold() must be positive"
        exit 105
    }

    if `p' < 1 {
        di as error "p() must be >= 1"
        exit 106
    }

    if `alpha' <= 0 | `alpha' >= 1 {
        di as error "alpha() must be in (0, 1)"
        exit 107
    }

    // Treatment binary check
    qui tab `treatment'
    if r(r) != 2 {
        di as error "treatment() must be binary (0/1)"
        exit 103
    }
    qui summarize `treatment'
    if r(min) != 0 | r(max) != 1 {
        di as error "treatment() must contain only 0 and 1"
        exit 103
    }

    // Return parsed results
    return local outcome "`varlist'"
    return local treatment "`treatment'"
    return local time "`time'"
    if "`cluster'" != "" {
        return local cluster "`cluster'"
    }
    return scalar threshold = `threshold'
    return scalar p = `p'
    return scalar alpha = `alpha'
    return scalar simulate = `simulate'
    return scalar seed = `seed'
    return scalar treat_time = `treat_time'
    return scalar nograph = ("`nograph'" != "")
    return scalar overall = ("`overall'" != "")
    return local mode = cond("`overall'" != "", "overall", "iterative")
end


