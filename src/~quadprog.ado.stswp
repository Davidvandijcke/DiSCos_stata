// quadprog.ado
program define quadprog, rclass
    version 14.0
    args mG mg0 mCE mce0 mCI mci0 mresults
    
    // Remove matrix prefix if present
    local G = subinstr("`mG'", "matrix(", "", .)
    local G = subinstr("`G'", ")", "", .)
    local g0 = subinstr("`mg0'", "matrix(", "", .)
    local g0 = subinstr("`g0'", ")", "", .)
    local CE = subinstr("`mCE'", "matrix(", "", .)
    local CE = subinstr("`CE'", ")", "", .)
    local ce0 = subinstr("`mce0'", "matrix(", "", .)
    local ce0 = subinstr("`ce0'", ")", "", .)
    local CI = subinstr("`mCI'", "matrix(", "", .)
    local CI = subinstr("`CI'", ")", "", .)
    local ci0 = subinstr("`mci0'", "matrix(", "", .)
    local ci0 = subinstr("`ci0'", ")", "", .)
    local results = subinstr("`mresults'", "matrix(", "", .)
    local results = subinstr("`results'", ")", "", .)
    
    plugin call quadprog_mata, `G' `g0' `CE' `ce0' `CI' `ci0' `results'
end

program quadprog_mata, plugin
