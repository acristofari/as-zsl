// ------------------------------------------------------------------------
// 
// This file is part of AS-ZSL, which is a solver for the following lasso
// problem with zero-sum constraint:
// 
//                     min 0.5*||Ax-y||^2 + lambda*||x||_1
//                s.t. sum(x) = 0
// 
// with given matrix A, vector y and scalar lambda.
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
// Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
// 
// Last update of this file:
// April 19th, 2022
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

#include "mex.h"
#include "as_zsl.h"
#include <string>

#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
typedef double mxDouble;
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    unsigned int m,n,n_lambda,t;
    std::vector<double> lambda;
    const double *cdbl_ptr;
    double *A,*y;
    mxDouble *mdbl_ptr;
    size_t *irs,*jcs;

    // check the number of inupts and outputs
    if (nrhs<3) {
        mexErrMsgTxt("At least three inputs are required.");
    }
    if (nrhs>4) {
        mexErrMsgTxt("At most four inputs are required.");
    }
    if (nlhs<1) {
        mexErrMsgTxt("At least one input is required.");
    }
    if (nlhs>3) {
        mexErrMsgTxt("At most outputs are required.");
    }

    // check inputs
    if (mxIsScalar(prhs[0]) || mxGetNumberOfDimensions(prhs[0])>2 || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("The first input must be a real matrix.");
    }
    m = (unsigned int) mxGetM(prhs[0]); // number of samples
    n = (unsigned int) mxGetN(prhs[0]); // number of variables
    if (mxGetNumberOfDimensions(prhs[0])>2 || mxGetN(prhs[1])!=1 || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgTxt("The second input must be a real column vector.");
    }
    if ((unsigned int) mxGetM(prhs[1]) != m) {
            mexErrMsgTxt("The dimension of A and y in the system Ax = y must agree.");
        }
    if (mxGetNumberOfDimensions(prhs[0])>2 || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgTxt("The third input must be a non-negative number or a column vector sorted in descending order.");
    }

    // get matrix A
    A = mxGetDoubles(prhs[0]);
    if (mxIsSparse(prhs[0])) {
        irs = mxGetIr(prhs[0]);
        jcs = mxGetJc(prhs[0]);
        for (unsigned int i=0; i<jcs[n]; i++) {
            if (mxIsNaN(A[i])) {
                mexErrMsgTxt("The first input must be a real matrix (NaNs are not allowed).");
            }
        }
    } else {
        for (unsigned int i=0; i<m*n; i++) {
            if (mxIsNaN(A[i])) {
                mexErrMsgTxt("The first input must be a real matrix (NaNs are not allowed).");
            }
        }
        irs = jcs = NULL;
    }

    // get vector y
    y = mxGetDoubles(prhs[1]);
    for (unsigned int i=0; i<m; i++) {
        if (mxIsNaN(y[i])) {
            mexErrMsgTxt("The second input must be a real column vector (NaNs are not allowed).");
        }
    }

    // get lambda(s)
    n_lambda = (unsigned int) mxGetN(prhs[2]);
    if (n_lambda > 1) {
        if (mxGetM(prhs[2]) > 1) {
            mexErrMsgTxt("The third input must be a non-negative number or a column vector sorted in descending order.");
        }
    } else {
        n_lambda = (unsigned int) mxGetM(prhs[2]);
    }
    lambda.resize(n_lambda);
    cdbl_ptr = mxGetDoubles(prhs[2]);
    lambda[0] = *cdbl_ptr;
    if (lambda[0]<0 || mxIsNaN(lambda[0])) {
        mexErrMsgTxt("The third input must be a non-negative number or a column vector sorted in descending order (NaNs are not allowed).");
    }
    for (unsigned int r=1; r<n_lambda; r++) {
        lambda[r] = cdbl_ptr[r];
        if (lambda[r]<0 || mxIsNaN(lambda[r]) || lambda[r]>lambda[r-1]) {
            mexErrMsgTxt("The third input must be a non-negative number or a column vector sorted in descending order (NaNs are not allowed).");
        }
    }
    
    // get options
    as_zsl_options opts;
    if (nrhs > 3) {
        if (!mxIsStruct(prhs[3]) || mxGetNumberOfElements(prhs[3])>1) {
            mexErrMsgTxt("The fourth input argument (which is optional) must be a structure.");
        }
        for (int i=0; i<mxGetNumberOfFields(prhs[3]); i++) {
            mxArray *tmp_mxArray = mxGetFieldByNumber(prhs[3],0,i);
            const char *tmp_char = mxGetFieldNameByNumber(prhs[3],i);
            if (std::string(tmp_char).compare(std::string("eps_opt")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'eps_opt' must be a non-negative number.");
                }
                opts.eps_opt = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("max_it")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'max_it' must be a number greater than or equal to 1.");
                }
                opts.max_it = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("verbosity")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsLogical(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'verbosity' must a logical.");
                }
                opts.verbosity = *mxGetLogicals(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("intercept")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsLogical(tmp_mxArray)) {
                    mexErrMsgTxt("In the options, 'intercept' must a logical.");
                }
                opts.intercept = *mxGetLogicals(tmp_mxArray);
            } else {
                mexErrMsgTxt("Not valid field name in the structure of options.");
            }
        }
    }

    // call the solver
    As_zsl alg(m,n,A,irs,jcs,y,lambda,&opts);
    alg.solve();
    
    // set outputs
    plhs[0] = mxCreateDoubleMatrix(n,n_lambda,mxREAL);
    mdbl_ptr = mxGetDoubles(plhs[0]);    
    t = 0;
    for (unsigned int r=0; r<n_lambda; r++) {
        cdbl_ptr = &alg.get_x()[r][0];
        for (unsigned int h=0; h<n; h++) {
            mdbl_ptr[t] = cdbl_ptr[h];
            t++;
        }
    }
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1,n_lambda,mxREAL);
        mdbl_ptr = mxGetDoubles(plhs[1]);
        cdbl_ptr = &alg.get_x0()[0];
        for (unsigned int r=0; r<n_lambda; r++) {
            mdbl_ptr[r] = cdbl_ptr[r];
        }
        if (nlhs > 2) {
            mxArray *pf = mxCreateDoubleMatrix(1,n_lambda,mxREAL);
            mxArray *pit = mxCreateDoubleMatrix(1,n_lambda,mxREAL);
            mxArray *pflag = mxCreateDoubleMatrix(1,n_lambda,mxREAL);
            const unsigned int *cuint_ptr;
            const char* as_zsl_field_names[] = {"f","it","flag"};
            mwSize as_zsl_info_dims[2] = {1,1};
            plhs[2] = mxCreateStructArray(2,as_zsl_info_dims,3,as_zsl_field_names);
            mdbl_ptr = mxGetDoubles(pf);
            cdbl_ptr = &alg.get_f()[0];
            for (unsigned int r=0; r<n_lambda; r++) {
                mdbl_ptr[r] = cdbl_ptr[r];
            }
            mxSetField(plhs[2],0,as_zsl_field_names[0],pf);
            mdbl_ptr = mxGetDoubles(pit);
            cuint_ptr = &alg.get_it()[0];
            for (unsigned int r=0; r<n_lambda; r++) {
                mdbl_ptr[r] = (double) cuint_ptr[r];
            }
            mxSetField(plhs[2],0,as_zsl_field_names[1],pit);
            mdbl_ptr = mxGetDoubles(pflag);
            cuint_ptr = &alg.get_flag()[0];
            for (unsigned int r=0; r<n_lambda; r++) {
                mdbl_ptr[r] = (double) cuint_ptr[r];
            }
            mxSetField(plhs[2],0,as_zsl_field_names[2],pflag);
        }
    }
    
}