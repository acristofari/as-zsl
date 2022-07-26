%  ------------------------------------------------------------------------
%  
%  This file is part of AS-ZSL, which is a solver for the following lasso
%  problem with zero-sum constraint:
%  
%                      min 0.5*||Ax-y||^2 + lambda*||x||_1
%                 s.t. sum(x) = 0
%  
%  with given matrix A, vector y and non-negative scalar lambda.
%  
%  ------------------------------------------------------------------------
%  
%  Reference paper:
%  
%  A. Cristofari (2022). A decomposition method for lasso problems with
%  zero-sum constraint. arXiv preprint 2204.07065
%  
%  ------------------------------------------------------------------------
%  
%  Author:
%  Andrea Cristofari (e-mail: andrea.cristofari@uniroma2.it)
%  
%  Last update of this file:
%  July 26th, 2022
% 
%  Licensing:
%  This file is part of AS-ZSL.
%  AS-ZSL is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  AS-ZSL is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.
%  You should have received a copy of the GNU General Public License
%  along with AS-ZSL. If not, see <http://www.gnu.org/licenses/>.
%  
%  Copyright 2022 Andrea Cristofari.
% 
%  ------------------------------------------------------------------------

clear all, clc;

rng(1);

% In this file, it is shown how to call AS-ZSL to solve a user-defined problem.

make; % build the MEX file (just the first time)

% (1) Get the problem (i.e., matrix of covariates 'A', response vector 'y' and regularization parameter 'lambda')
m = 1000; % number of samples
n = 5000; % number of features
E = exp(randn(m,n));
A = log(E./sum(E,2)); % make data compositional and then apply a log transformation
% A = A - mean(A); % to center the columns of A
y = A*sprand(n,1,0.05) + normrnd(0,0.5,m,1);
% y = y - mean(y); % to center y
lambda = 0.5*peak2peak(A'*y)/2;

% (2) Call AS-ZSL
[x,f,as_zsl_info] = as_zsl(A,y,lambda);

% If 'lambda' is a vector of regularization parameters sorted in descending order,
% then we have to solve a problem for each regularization parameter. In this case,
% a warm start strategy is used and the output values have columns, each one
% referring to the corresponding regularization parameter.
% For example:
% 
% lambda_max = peak2peak(A'*y)/2;
% lambda = logspace(log10(0.95*lambda_max),log10(1e-3*lambda_max),5);
% [x,f,as_zsl_info] = as_zsl(A,y,lambda);


%--------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE AS-ZSL PARAMETERS ***
% (see the file 'usage.txt' to know which parameters can be changed and
% their default values)
%
% Instead of calling AS-ZSL by the above instruction, do the following:
%
% - create a structure having as field names the names of the parameters
%   to be changed and assign them new values, e.g.,
%
%     opts.verbosity = true;
%
% - pass the structure to AS-ZSL as fourth input argument, e.g.,
%
%     [x,f,as_zsl_info] = as_zsl(A,y,lambda,opts);
%--------------------------------------------------------------------------