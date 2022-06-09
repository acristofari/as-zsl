#ifndef __AS_ZSL_H_INCLUDED__
#define __AS_ZSL_H_INCLUDED__
#include <vector>
#include <cstddef>

struct as_zsl_options {

    // *** DO NOT CHANGE DEFINITIONS IN THIS STRUCTURE ***

    // ==============================================
    // DEFAULT VALUES OF AS-ZSL PARAMETERS
    // ==============================================

    bool intercept = true;
    double eps_opt = 1e-3;
    int max_it = 100000;
    bool verbosity = false;

};

class As_zsl {
private:
    bool intercept,A_is_full;
    unsigned int m,n,n_lambda,it,max_it,verbosity,n_non_act,flag,i,j;
    double lambda,eps_opt,f,f_quad,f_norm1,g_i,g_j;
    size_t *irs,*jcs;
    const double *lambda_vec;
    std::vector<bool> is_fixed;
    std::vector<unsigned int> it_vec,flag_vec,ind_non_act;
    std::vector<double> x,u,f_vec,x0_vec;
    std::vector<std::vector<double>> x_vec;
    double *A,*y;
    void solve_subproblem();

public:
    As_zsl(unsigned int, unsigned int, double * const, size_t * const, size_t * const,
        double * const, const std::vector<double>&, const as_zsl_options*);
    As_zsl(unsigned int n_row, unsigned int n_col, double * const mat, size_t * const row, size_t * const col,
        double * const label, const std::vector<double>& lam) : As_zsl(n_row, n_row, mat, row, col, label, lam, NULL){};
    void solve();
    const std::vector<std::vector<double>>& get_x();
    const std::vector<double>& get_x0();
    const std::vector<double>& get_f();
    const std::vector<unsigned int>& get_it();
    const std::vector<unsigned int>& get_flag();
};

#endif