program define disco_weight, rclass
    version 18.0
    
    // Check if disco was run
    if "`e(cmd)'" != "disco" {
        di as error "disco_weight is a post-estimation command for disco. Run disco first."
        exit 198
    }
    
    // Syntax: disco_weight id_var name_var [, n(#) format(str) frame(str) round(#)]
    syntax varlist(min=2 max=2) [, N(integer 5) Format(string) FRAMe(string) Round(real 0.0001)]
    
    // Parse variables
    tokenize `varlist'
    local id_var "`1'"
    local name_var "`2'"
    
    // Validate n
    if `n' <= 0 {
        di as error "n() must be positive"
        exit 198
    }
    
    // Default format
    if "`format'" == "" local format "%12.4f"
    
    // Default frame name if not specified
    if "`frame'" == "" local frame "disco_weights"
    
    // Create a new frame for results
    cap frame drop `frame'
    frame create `frame' str32 name double weight
    frame `frame' {
        // Loop through cids and weights
        forvalues i = 1/`=colsof(e(cids))' {
            local id = e(cids)[1,`i']
            local w = e(weights)[1,`i']
            qui frame default: levelsof `name_var' if `id_var' == `id', local(name) clean
            qui set obs `=_N + 1'
            qui replace name = "`name'" if _n == _N
            qui replace weight = `w' if _n == _N
        }
        
        // Sort and keep top n
        gsort -weight
        qui keep in 1/`n'
        replace weight = round(weight, `round')
    }
    
    // Display results
    di _n as txt "Top `n' weights:"
    di as txt "{hline 50}"
    frame `frame': list name weight, noobs table div
    
    // Return results
    return clear
    return local cmd "disco_weight"
    return scalar n = `n'
    return local format "`format'"
    return local frame "`frame'"
end
