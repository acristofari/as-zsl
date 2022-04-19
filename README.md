# A decomposition method for lasso problems with zero-sum constraint

_Active-Set Zero-Sum-Lasso_ (AS-ZSL) is a solver for the following lasso problem with zero-sum constraint:

         min 0.5*||Ax-y||^2 + lambda*||x||_1
    s.t. sum(x) = 0

with given matrix _A_, vector _y_ and scalar _lambda_.

AS-ZSL combines a tailored _active-set technique_, to identify the zero variables in the optimal solution,
with a _2-coordinate descent scheme_.

This software is written in C++ and can be called from Matlab using a MEX file.

## Reference paper

[A. Cristofari (2022). _A decomposition method for lasso problems with zero-sum constraint._ arXiv preprint 2204.07065](https://arxiv.org/abs/2204.07065).


## Author

Andrea Cristofari (e-mail: [andrea.cristofari@unipd.it](mailto:andrea.cristofari@unipd.it))

## Licensing

AS-ZSL is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
AS-ZSL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with AS-ZSL. If not, see <http://www.gnu.org/licenses/>.

Copyright 2022 Andrea Cristofari.

## How to use AS-ZSL in Matlab

1. This directory should contain the following files:

    * `as_zsl.cpp`,
    * `as_zsl.h`,
    * `as_zsl_matlab.cpp`,
    * `COPYING.txt`,
    * `main.m`,
    * `make.m`,
    * `README.md`,
    * `usage.txt`.

2. In Matlab, run `make.m` to build the MEX file.

3. See the file `usage.txt` to know how to call AS-ZSL from Matlab, change algorithm parameters and get output values.

4. See the file `main.m` for an example. To run the example, just call `main.m` in Matlab.
