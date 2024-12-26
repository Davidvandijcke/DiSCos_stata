// clear
// matrix G = (2,0 \ 0,2)
// matrix g0 = (-2 \ -5)
// matrix CE = (1 \ 1)
// matrix ce0 = (1)
// matrix CI = (1,0 \ 0,1)
// matrix ci0 = (0 \ 0)
// matrix results = J(1, 3, 0)  // room for 2 variables + cost
//
// gen y = 1
// plugin call quadprog_mata, G g0 CE ce0 CI ci0 results
//

// Load everything

program drop _all

mata: mata clear
clear

do quadprog.ado

do quadprog.mata


// Define matrices
matrix G = (1,0 \ 0,1) // Identity matrix
matrix g0 = (0 \ 0)    // 2x1 zero vector
matrix CE = (1 \ 1)    // 2x1
matrix ce0 = (-1)      // p=1, so ce0 is 1x1 with value -1
matrix CI = (1,0 \ 0,1) // This enforces x1 >=0 and x2 >=0
matrix ci0 = (0 \ 0)    // 2x1 zero vector
matrix results = J(1, 3, 0)  // 1 row, 3 columns for x1, x2, cost

// Create a dummy variable to ensure the plugin call knows there's data
gen y = 1

// // Run the quadprog interface defined in quadprog.ado and quadprog.mata
// mata:
// G_m = st_matrix("G")
// g0_m = st_matrix("g0")
// CE_m = st_matrix("CE")
// ce0_m = st_matrix("ce0")
// CI_m = st_matrix("CI")
// ci0_m = st_matrix("ci0")
//
// res = solve_quadprog(G_m, g0_m, CE_m, ce0_m, CI_m, ci0_m)
// res
// end

// The above 'res' should contain a row vector: [x1, x2, cost] = [0.5, 0.5, 0.25]

// Or call via quadprog ado directly
quadprog matrix(G) matrix(g0) matrix(CE) matrix(ce0) matrix(CI) matrix(ci0) matrix(results)

// Check the results matrix
matrix list results
