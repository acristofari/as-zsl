// ------------------------------------------------------------------------
// 
// This file is part of AS-ZSL, which is a solver for the following lasso
// problem with zero-sum constraint:
// 
//                     min 0.5*||Ax-y||^2 + lambda*||x||_1
//                s.t. sum(x) = 0
// 
// with given matrix A, vector y and non-negative scalar lambda.
// 
// ------------------------------------------------------------------------
// 
// Reference paper:
// 
// A. Cristofari (2022). A decomposition method for lasso problems with
// zero-sum constraint. arXiv preprint 2204.07065
// 
// ------------------------------------------------------------------------
// 
// Author:
// Andrea Cristofari (e-mail: andrea.cristofari@uniroma2.it)
// 
// Last update of this file:
// July 26th, 2022
//
// Licensing:
// This file is part of AS-ZSL.
// AS-ZSL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// AS-ZSL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with AS-ZSL. If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2022 Andrea Cristofari.
//
// ------------------------------------------------------------------------

#include "as_zsl.h"
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>

// constructor
//-------------------------------------------------------------------------------------
As_zsl::As_zsl(unsigned int n_row, unsigned int n_col, double * const mat,
    size_t * const row, size_t * const col, double * const label,  const std::vector<double>& lam,
    const as_zsl_options * opts) {

    if (opts == NULL) {
        as_zsl_options as_zsl_opts;
        opts = &as_zsl_opts;
    }

    eps_opt = opts->eps_opt;
    max_it = (unsigned int) opts->max_it;
    verbosity = opts->verbosity;

    if (eps_opt < 0e0) {
        std::cout << "In the options, 'eps_opt' must be a non-negative number.\n";
        exit(-1);
    }
    if (opts->max_it < 1) {
        std::cout << "In the options, 'max_it' must be a number greater than or equal to 1.\n";
        exit(-1);
    }

    m = n_row; // number of columns of A
    n = n_col; // number of rows of A
    A = mat;
    y = label;
    lambda_vec = &lam[0];
    n_lambda = (unsigned int) lam.size();

    x_vec.resize(n_lambda);
    f_vec.resize(n_lambda);
    it_vec.resize(n_lambda);
    flag_vec.resize(n_lambda);

    if (row==NULL || col==NULL) {
        if (row!=NULL || col!=NULL) {
            std::cout << "The matrix of covariates must be either full or sparse.\n";
            exit(-1);
        }
        A_is_full = true;
    } else {
        A_is_full = false;
        irs = row;
        jcs = col;
    }
    
}
//-------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------
void As_zsl::solve() {

    unsigned int seed,h_non_act;
    double opt_viol,eta_min,eta_max,f_old,mu,theta,tmp;
    double *A_col;
    bool it_mvp;
    std::vector<double> g,c,pi;
    
    if (verbosity) {
        std::cout.precision(4);
        std::cout.setf(std::ios::scientific,std::ios::floatfield);
    }
    
    is_fixed.assign(n,false);
    c.resize(n);
    g.resize(n);
    pi.resize(n);
    // set c = A'*y
    if (A_is_full) {
        for (unsigned int h=0; h<n; h++) {
            tmp = 0e0;
            A_col = A + h*m;
            for (unsigned int t=0; t<m; t++) {
                tmp += A_col[t]*y[t];
            }
            c[h] = tmp;
        }
    } else {
        for (unsigned int h=0; h<n; h++) {
            tmp = 0e0;
            for (size_t t=jcs[h]; t<jcs[h+1]; t++) {
                tmp += A[t]*y[irs[t]];
            }
            c[h] = tmp;
        }
    }

    // non-active indices
    seed = 1;
    srand(seed);
    std::minstd_rand0 eng(seed);
    n_non_act = n;
    ind_non_act.resize(n);
    std::iota(ind_non_act.begin(),ind_non_act.end(),0);

    for (unsigned int r=0; r<n_lambda; r++) {

        x_vec[r].resize(n);
        lambda = lambda_vec[r];

        if (verbosity && n_lambda>1) {
            std::cout << "\n*** lambda = " << lambda << " ***\n";
        }

        // warm start in case of multiple values of lambda
        if (r > 0) {

            if (f_norm1 > 0e0) {

                f = f_quad + lambda*f_norm1;
                it = 0;
                if (verbosity) {
                    std::cout << "it = 0, f = " << f << "\n";
                }

            } else {

                if (verbosity) {
                    std::cout << "it = 0, f = " << f << "\n";
                }

                eta_min -= lambda_vec[r-1];
                eta_min += lambda;
                eta_max += lambda_vec[r-1];
                eta_max -= lambda;
                opt_viol =  eta_max - eta_min;
                if (opt_viol <= eps_opt) {
                    x_vec[r] = x;
                    f_vec[r] = f;
                    it_vec[r] = 0;
                    flag_vec[r] = 0;
                    continue;
                }
                
                // first mvp iteration
                solve_subproblem();
                if (verbosity) {
                    std::cout << "it = 1, f = " << f << " (mvp)\n";
                }

                it = 1;

            }

        } else {
            
            x.assign(n,0e0);
            f_quad = 0e0;
            for (unsigned int h=0; h<m; h++) {
                f_quad += y[h]*y[h];
            }
            f_quad *= 5e-1;
            f = f_quad;
            f_norm1 = 0e0;
            for (unsigned int h=0; h<n; h++) {
                g[h] = -c[h];
            }
            u.assign(m,0e0);

            if (verbosity) {
                std::cout << "it = 0, f = " << f << "\n";
            }

            // compute the mvp
            i = j = 0;
            eta_min = eta_max = g[0];
            for (unsigned int h=1; h<n; h++) {
                if (g[h] < eta_min) {
                    eta_min = g[h];
                    i = h;
                } else if (g[h] > eta_max) {
                    eta_max = g[h];
                    j = h;
                }
            }
            eta_min += lambda;
            eta_max -= lambda;
            opt_viol =  eta_max - eta_min;
            if (opt_viol <= eps_opt) {
                x_vec[r] = x;
                f_vec[r] = f;
                it_vec[r] = 0;
                flag_vec[r] = 0;
                continue;
            }
            
            // first mvp iteration
            g_i = g[i];
            g_j = g[j];
            solve_subproblem();
            if (verbosity) {
                std::cout << "it = 1, f = " << f << " (mvp)\n";
            }

            it = 1;

        }

        theta = 1e-2;
        flag = 1;

        // fake values for the 1st iteration
        it_mvp = false;
        f_old = 0e0;

        while (it < max_it) {
            
            // choose between strategy mvp and strategy ac2cd
            if (!it_mvp && f_old-f<=theta*std::max(f_old,1e0)) {
                it_mvp = true;
                theta = std::max(5e-1*theta,1e-6);
            } else {
                if (it_mvp && f_old-f<=1e-9*std::max(f_old,1e0)) {
                    flag = 2;
                    break;
                }
                it_mvp = false;
            }
        
            if (it_mvp) {

                // compute the gradient
                if (A_is_full) {
                    for (unsigned int h=0; h<n; h++) {
                        tmp = 0e0;
                        A_col = A + h*m;
                        for (unsigned int t=0; t<m; t++) {
                            tmp += A_col[t]*u[t];
                        }
                        g[h] = tmp - c[h];
                    }
                } else {
                    for (unsigned int h=0; h<n; h++) {
                        tmp = 0e0;
                        for (size_t t=jcs[h]; t<jcs[h+1]; t++) {
                            tmp += A[t]*u[irs[t]];
                        }
                        g[h] = tmp - c[h];
                    }
                }

                // compute the multiplier
                mu = 0e0;
                tmp = 0e0;
                for (unsigned int h=0; h<n; h++) {
                    if (x[h]!=0e0) {
                        mu += fabs(x[h])*(g[h]+(x[h]>0e0?lambda:-lambda));
                        tmp += fabs(x[h]);
                    }
                }
                mu /= tmp;
                for (unsigned int h=0; h<n; h++) {
                    pi[h] = g[h] - mu;
                }

            }
        
            // active-set estimate
            n_non_act = 0;
            for (unsigned int h=0; h<n; h++) {
                if (x[h]!=0 || (fabs(pi[h])>lambda && !is_fixed[i])) {
                    ind_non_act[n_non_act] = h;
                    n_non_act++;
                }
            }

            f_old = f;
        
            if (it_mvp) {

                eta_min = std::numeric_limits<double>::max();
                eta_max = std::numeric_limits<double>::lowest();

                for (unsigned int h=0; h<n_non_act; h++) {
                    h_non_act = ind_non_act[h];
                    if (x[h_non_act] == 0e0) {
                        if (g[h_non_act]+lambda < eta_min) {
                            eta_min = g[h_non_act] + lambda;
                            i = h_non_act;
                        }
                        if (g[h_non_act]-lambda > eta_max) {
                            eta_max = g[h_non_act] - lambda;
                            j = h_non_act;
                        }
                    } else {
                        tmp = g[h_non_act] + (x[h_non_act]>0e0?lambda:-lambda);
                        if (tmp < eta_min) {
                            eta_min = tmp;
                            i = h_non_act;
                        }
                        if (tmp > eta_max) {
                            eta_max = tmp;
                            j = h_non_act;
                        }
                    }
                }

                opt_viol =  eta_max - eta_min; // this is the optimality violation in N^k

                if (opt_viol <= eps_opt) {
                    // compute the overall optimality violation
                    if (x[0] == 0e0) {
                        eta_min = g[0] + lambda;
                        eta_max = g[0] - lambda;
                    } else {
                        eta_min = eta_max = g[0] + (x[0]>0?lambda:-lambda);
                    }
                    for (unsigned int h=1; h<n; h++) {
                        if (x[h] == 0e0) {
                            if (g[h]+lambda < eta_min) {
                                eta_min = g[h] + lambda;
                                i = h;
                            }
                            if (g[h]-lambda > eta_max) {
                                eta_max = g[h] - lambda;
                                j = h;
                            }
                        } else {
                            tmp = g[h] + (x[h]>0e0?lambda:-lambda);
                            if (tmp < eta_min) {
                                eta_min = tmp;
                                i = h;
                            }
                            if (tmp > eta_max) {
                                eta_max = tmp;
                                j = h;
                            }
                        }
                    }
                }

                opt_viol = eta_max - eta_min;
                if (opt_viol <= eps_opt) {
                    flag = 0;
                    break;
                }
                        
                g_i = g[i];
                g_j = g[j];

                solve_subproblem();

            } else {

                // shuffle variables
                shuffle(ind_non_act.begin(),ind_non_act.begin()+n_non_act,std::default_random_engine(seed=eng()));

                // select index j
                tmp = fabs(x[ind_non_act[0]]);
                j = 0;
                for (unsigned int h=1; h<n_non_act; h++) {
                    if (fabs(x[ind_non_act[h]]) > tmp) {
                        tmp = fabs(x[ind_non_act[h]]);
                        j = h;
                    }
                }
            
                // solve subproblems
                j = ind_non_act[j];
                ind_non_act[j] = ind_non_act[0];
                for (unsigned int h=1; h<n_non_act; h++) {
                    i = ind_non_act[h];
                    if (A_is_full) {
                        tmp = 0e0;
                        A_col = A + i*m;
                        for (unsigned int t=0; t<m; t++) {
                            tmp += A_col[t]*u[t];
                        }
                        g_i = tmp - c[i];
                        tmp = 0e0;
                        A_col = A + j*m;
                        for (unsigned int t=0; t<m; t++) {
                            tmp += A_col[t]*u[t];
                        }
                        g_j = tmp - c[j];
                    } else {
                        tmp = 0e0;
                        for (size_t t=jcs[i]; t<jcs[i+1]; t++) {
                            tmp += A[t]*u[irs[t]];
                        }
                        g_i = tmp - c[i];
                        tmp = 0e0;
                        for (size_t t=jcs[j]; t<jcs[j+1]; t++) {
                            tmp += A[t]*u[irs[t]];
                        }
                        g_j = tmp - c[j];
                    }
                    solve_subproblem();
                }
            
            }

            it++;

            if (verbosity) {
                std::cout << "it = " << it << ", f = " << f;
                if (it_mvp) {
                    std::cout << " (mvp)\n";
                } else {
                    std::cout << " (ac2cd)\n";
                }
            }

        }

        x_vec[r] = x;
        f_vec[r] = f;
        it_vec[r] = it;
        flag_vec[r] = flag;

    }

}
//-------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------
void As_zsl::solve_subproblem() {

    double x_i,alpha,beta,s;
    std::vector<double> A_ij;
    double sum_xi_xj = x[i] + x[j];
    bool sol_found = false;

    alpha = 0e0;
    if (A_is_full) {
        A_ij.resize(n);
        double *A_col_i = A + m*i;
        double *A_col_j = A + m*j;
        for (unsigned int h=0; h<m; h++) {
            A_ij[h] = A_col_i[h] - A_col_j[h];
            alpha += A_ij[h]*A_ij[h];
        }
    } else {
        A_ij.assign(n,0e0);
        for (size_t h=jcs[i]; h<jcs[i+1]; h++) {
            A_ij[irs[h]] = A[h];
        }
        for (size_t h=jcs[j]; h<jcs[j+1]; h++) {
            A_ij[irs[h]] -= A[h];
        }
        alpha = std::inner_product(A_ij.begin(),A_ij.end(),A_ij.begin(),0e0);
    }
    
    if (alpha > 0e0) {

        beta = alpha*x[i] + g_j - g_i;

        // seek a stationary point
        x_i = (beta-2e0*lambda)/alpha;
        if (x_i > std::max(0e0,sum_xi_xj)) {
            sol_found = true;
        } else {
            x_i = (beta+2e0*lambda)/alpha;
            if (x_i < std::min(0e0,sum_xi_xj)) {
                sol_found = true;
            } else {
                x_i = beta/alpha;
                if (x_i*(x_i-sum_xi_xj)<0e0) {
                    sol_found = true;
                }
            }
        }

        // stationary point not found -> the minimizer is a point of non-differentiability
        if (!sol_found) {
            if (lambda*fabs(sum_xi_xj) <= 5e-1*alpha*sum_xi_xj*sum_xi_xj-beta*sum_xi_xj+lambda*fabs(sum_xi_xj)) {
                x_i = 0e0;
            } else {
                x_i = sum_xi_xj;
            }
        }

    } else {
        
        x_i = 0e0;
        is_fixed[i] = true;
        
    }

    // update x, f and u
    s = x_i - x[i];
    if (s != 0e0) {
        f_quad += s*(g_i-g_j+5e-1*s*alpha);
        f_norm1 += fabs(x_i) + fabs(x[j]-s) - fabs(x[i]) - fabs(x[j]);
        f = f_quad + lambda*f_norm1;
        for (unsigned int h=0; h<m; h++) {
            u[h] += s*A_ij[h];
        }
        x[i] = x_i;
        x[j] = x[j] - s;
    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<std::vector<double>>& As_zsl::get_x() {
    return x_vec;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<double>& As_zsl::get_f() {
    return f_vec;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<unsigned int>& As_zsl::get_it() {
    return it_vec;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<unsigned int>& As_zsl::get_flag() {
    return flag_vec;
}
//-------------------------------------------------------------------------------------