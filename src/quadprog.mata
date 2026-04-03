// quadprog.mata
version 18.0
mata:

matrix solve_quadprog(real matrix G, 
                                    real vector g0, 
                                    real matrix CE, 
                                    real vector ce0, 
                                    real matrix CI, 
                                    real vector ci0)
{
    real scalar n
    string scalar cmd
    string scalar tmpG, tmpg0, tmpCE, tmpce0, tmpCI, tmpci0, tmpres
    real vector res
    
    // Get dimensions
    n = rows(G)
    
    // Create temporary names for Stata matrices
    tmpG = st_tempname()
    tmpg0 = st_tempname()
    tmpCE = st_tempname()
    tmpce0 = st_tempname()
    tmpCI = st_tempname()
    tmpci0 = st_tempname()
    tmpres = st_tempname()
    
    // Transfer matrices to Stata
    st_matrix(tmpG, G)
    st_matrix(tmpg0, g0)
    st_matrix(tmpCE, CE)
    st_matrix(tmpce0, ce0)
    st_matrix(tmpCI, CI)
    st_matrix(tmpci0, ci0)
    
    // Initialize results matrix
    st_matrix(tmpres, J(1, n + 1, 0))
    
    // Build command string
    cmd = "quadprog " + 
          "matrix(" + tmpG + ") " + 
          "matrix(" + tmpg0 + ") " + 
          "matrix(" + tmpCE + ") " + 
          "matrix(" + tmpce0 + ") " + 
          "matrix(" + tmpCI + ") " + 
          "matrix(" + tmpci0 + ") " + 
          "matrix(" + tmpres + ")"
    
    // Print command for debugging
    // printf("Executing command: %s\n", cmd)
    
    // Call Stata command
    stata(cmd)
    
    // Get results back
    res = st_matrix(tmpres)
    
    
    return(res)
}

end
