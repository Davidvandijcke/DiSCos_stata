// quadprog_mata.cpp
#include <iostream>
#include <vector>
#include "Array.hh"      // Must come before QuadProg++.hh
#include "QuadProg++.hh"
#include "stplugin.h"
#include <string>

using namespace quadprogpp;

ST_plugin *_stata_; // define the global _stata_ variable
#define ERROR_MSG(x) (const_cast<char*>(x))

STDLL pginit(ST_plugin *p)
{
    _stata_ = p;
    return(SD_PLUGINVER);
}

STDLL stata_call(int argc, char *argv[])
{
    // We expect 7 arguments: G g0 CE ce0 CI ci0 results
    if (argc < 7) {
        std::string error_message = "Not enough arguments provided.";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }

    // argv[0] = G
    // argv[1] = g0
    // argv[2] = CE
    // argv[3] = ce0
    // argv[4] = CI
    // argv[5] = ci0
    // argv[6] = results

    // Get dimensions of G
    ST_int n_vars = SF_row(argv[0]);
    ST_int n_vars2 = SF_col(argv[0]);
    if (n_vars == 0 || n_vars2 == 0) {
        std::string error_message = "Matrix G not found or invalid.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }
    if (n_vars != n_vars2) {
        std::string error_message = "G must be a square matrix.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }

    // g0 must be n_vars x 1
    ST_int g0_rows = SF_row(argv[1]);
    ST_int g0_cols = SF_col(argv[1]);
    if (g0_rows != n_vars || g0_cols != 1) {
        std::string error_message = "g0 must be an n x 1 vector.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }

    // CE must be n_vars x p
    ST_int CE_rows = SF_row(argv[2]);
    ST_int CE_cols = SF_col(argv[2]);
    if (CE_rows != n_vars) {
        std::string error_message = "CE must have same number of rows as G.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }

    // ce0 must be p x 1
    ST_int ce0_rows = SF_row(argv[3]);
    ST_int ce0_cols = SF_col(argv[3]);
    if (ce0_rows != CE_cols || ce0_cols != 1) {
        std::string error_message = "ce0 must be p x 1, where p = number of equality constraints.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }
    ST_int n_eq_constraints = CE_cols;

    // CI must be n_vars x m
    ST_int CI_rows = SF_row(argv[4]);
    ST_int CI_cols = SF_col(argv[4]);
    if (CI_rows != n_vars) {
        std::string error_message = "CI must have same number of rows as G.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }

    // ci0 must be m x 1
    ST_int ci0_rows = SF_row(argv[5]);
    ST_int ci0_cols = SF_col(argv[5]);
    if (ci0_rows != CI_cols || ci0_cols != 1) {
        std::string error_message = "ci0 must be m x 1, where m = number of inequality constraints.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }
    ST_int n_ineq_constraints = CI_cols;

    // results must have at least 1 row and n_vars+1 columns
    ST_int res_rows = SF_row(argv[6]);
    ST_int res_cols = SF_col(argv[6]);
    if (res_rows < 1 || res_cols < n_vars + 1) {
        std::string error_message = "results matrix must have at least 1 row and (n_vars+1) columns.\n";
        SF_error(const_cast<char*>(error_message.c_str()));
        return(198);
    }

    try {
        // Create temporary arrays
        std::vector<double> G_data(n_vars * n_vars);
        std::vector<double> g0_data(n_vars);
        std::vector<double> CE_data(n_vars * n_eq_constraints);
        std::vector<double> ce0_data(n_eq_constraints);
        std::vector<double> CI_data(n_vars * n_ineq_constraints);
        std::vector<double> ci0_data(n_ineq_constraints);

        double val;

        // Read G
        for (int i = 1; i <= n_vars; i++) {
            for (int j = 1; j <= n_vars; j++) {
                if (SF_mat_el(argv[0], i, j, &val)) {
                    std::string error_message = "Error reading G.\n";
                    SF_error(const_cast<char*>(error_message.c_str()));
                    return(198);
                }
                G_data[(i-1)*n_vars + (j-1)] = val;
            }
        }

        // Read g0
        for (int i = 1; i <= n_vars; i++) {
            if (SF_mat_el(argv[1], i, 1, &val)) {
                std::string error_message = "Error reading g0.\n";
                SF_error(const_cast<char*>(error_message.c_str()));
                return(198);
            }
            g0_data[i-1] = val;
        }

        // Read CE
        for (int i = 1; i <= n_vars; i++) {
            for (int j = 1; j <= n_eq_constraints; j++) {
                if (SF_mat_el(argv[2], i, j, &val)) {
                    std::string error_message = "Error reading CE.\n";
                    SF_error(const_cast<char*>(error_message.c_str()));
                    return(198);
                }
                CE_data[(i-1)*n_eq_constraints + (j-1)] = val;
            }
        }

        // Read ce0
        for (int i = 1; i <= n_eq_constraints; i++) {
            if (SF_mat_el(argv[3], i, 1, &val)) {
                std::string error_message = "Error reading ce0.\n";
                SF_error(const_cast<char*>(error_message.c_str()));
                return(198);
            }
            ce0_data[i-1] = val;
        }

        // Read CI
        for (int i = 1; i <= n_vars; i++) {
            for (int j = 1; j <= n_ineq_constraints; j++) {
                if (SF_mat_el(argv[4], i, j, &val)) {
                    std::string error_message = "Error reading CI.\n";
                    SF_error(const_cast<char*>(error_message.c_str()));
                    return(198);
                }
                CI_data[(i-1)*n_ineq_constraints + (j-1)] = val;
            }
        }

        // Read ci0
        for (int i = 1; i <= n_ineq_constraints; i++) {
            if (SF_mat_el(argv[5], i, 1, &val)) {
                std::string error_message = "Error reading ci0.\n";
                SF_error(const_cast<char*>(error_message.c_str()));
                return(198);
            }
            ci0_data[i-1] = val;
        }

        // Build QuadProg++ structures
        quadprogpp::Matrix<double> G(n_vars, n_vars);
        quadprogpp::Vector<double> g0(n_vars);
        quadprogpp::Matrix<double> CE(n_vars, n_eq_constraints);
        quadprogpp::Vector<double> ce0(n_eq_constraints);
        quadprogpp::Matrix<double> CI(n_vars, n_ineq_constraints);
        quadprogpp::Vector<double> ci0(n_ineq_constraints);
        quadprogpp::Vector<double> x(n_vars);

        for (int i = 0; i < n_vars; i++) {
            for (int j = 0; j < n_vars; j++) {
                G[i][j] = G_data[i*n_vars + j];
            }
            g0[i] = g0_data[i];
        }

        for (int i = 0; i < n_vars; i++) {
            for (int j = 0; j < n_eq_constraints; j++) {
                CE[i][j] = CE_data[i*n_eq_constraints + j];
            }
        }

        for (int i = 0; i < n_eq_constraints; i++) {
            ce0[i] = ce0_data[i];
        }

        for (int i = 0; i < n_vars; i++) {
            for (int j = 0; j < n_ineq_constraints; j++) {
                CI[i][j] = CI_data[i*n_ineq_constraints + j];
            }
        }

        for (int i = 0; i < n_ineq_constraints; i++) {
            ci0[i] = ci0_data[i];
        }

        // Solve QP
        double result = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

        // Store results: x and then cost
        for (int i = 0; i < n_vars; i++) {
            if (SF_mat_store(argv[6], 1, i+1, x[i])) {
                std::string error_message = "Error storing solution in results.\n";
                SF_error(const_cast<char*>(error_message.c_str()));
                return(198);
            }
        }

        if (SF_mat_store(argv[6], 1, n_vars+1, result)) {
            std::string error_message = "Error storing cost in results.\n";
            SF_error(const_cast<char*>(error_message.c_str()));
            return(198);
        }

        return(0);
    }
    catch (const std::exception &e) {
        SF_error(ERROR_MSG(e.what()));
        return(198);
    }
}

