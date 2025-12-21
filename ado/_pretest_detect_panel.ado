*! _pretest_detect_panel v0.1.0 - Panel Data Detection Module
*!
*! Description:
*!   Detects whether data is set as panel data via xtset.
*!   - If panel data, returns panelvar and timevar
*!   - If not panel data, marks as repeated cross-section
*!
*! Syntax:
*!   _pretest_detect_panel
*!
*! Returns:
*!   r(is_panel)  - 1 if panel data, 0 if repeated cross-section
*!   r(panelvar)  - Panel variable name (only if is_panel=1)
*!   r(timevar)   - Time variable name (only if is_panel=1)
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), Section 2.1:
*!   "Our results apply both to a panel data setting, where n subjects 
*!    have been tracked over the T time periods, or to repeated cross-sections"

program define _pretest_detect_panel, rclass
    version 17.0
    
    // Attempt to get xtset information
    capture xtset
    local xtset_rc = _rc
    
    if `xtset_rc' == 0 {
        // xtset is configured, this is panel data
        local panelvar = r(panelvar)
        local timevar = r(timevar)
        
        // Verify returned variable names are non-empty
        if "`panelvar'" != "" & "`panelvar'" != "." {
            return scalar is_panel = 1
            return local panelvar "`panelvar'"
            return local timevar "`timevar'"
            
            // Count unique individuals (n in paper)
            tempvar tag
            qui egen `tag' = tag(`panelvar')
            qui count if `tag' == 1
            local n_units = r(N)
            return scalar n_units = `n_units'
            
            // Compute balanced panel information
            // Note: Use levelsof to get actual period count, avoiding time gaps
            qui levelsof `timevar', local(tlevels)
            local n_periods : word count `tlevels'
            qui count
            local n_obs = r(N)
            local is_balanced = (`n_obs' == `n_units' * `n_periods')
            return scalar is_balanced = `is_balanced'
            return scalar n_periods = `n_periods'
        }
        else {
            // xtset returned empty panelvar (maybe only time was set)
            // Treat as repeated cross-section
            return scalar is_panel = 0
            return local panelvar ""
            return local timevar ""
        }
    }
    else {
        // xtset not configured, this is repeated cross-section data
        return scalar is_panel = 0
        return local panelvar ""
        return local timevar ""
    }
end
