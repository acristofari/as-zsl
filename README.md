# A decomposition method for lasso problems with zero-sum constraint

_Active-Set Zero-Sum-Lasso_ (AS-ZSL) is a solver for the following lasso problem with zero-sum constraint:

<img src="https://latex.codecogs.com/svg.image?\min&space;\frac12&space;||Ax-y||^2&space;&plus;&space;\lambda||x||_1&space;\\\text{s.t.&space;}&space;&space;\sum_{i=1}^n&space;x_i&space;=&space;0">

with given matrix A, vector y and non-negative scalar &lambda;.

AS-ZSL combines a tailored _active-set technique_, to identify the zero variables in the optimal solution,
with a _2-coordinate descent scheme_.

This software is written in C++ and can be called from Matlab using a MEX file.

## Reference paper

[A. Cristofari (2024). A decomposition method for lasso problems with zero-sum constraint. European Journal of Operational Research 306(1), 358–369](https://doi.org/10.1016/j.ejor.2022.09.030).

Author's note: due to a mistake during the publication process, problem (1) is erroneously referred to as (A.3) throughout the paper (an arXiv version is available which might help)

## Author

Andrea Cristofari (e-mail: [andrea.cristofari@uniroma2.it](mailto:andrea.cristofari@uniroma2.it))

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

Copyright 2022-2024 Andrea Cristofari.

## How to use AS-ZSL in Matlab

1. In Matlab, run `make.m` to build the MEX file.

2. See the file `usage.txt` to know how to call AS-ZSL from Matlab, change algorithm parameters and get output values.

3. See the file `main.m` for an example. To run the example, just call `main.m` in Matlab.
