# qmsm

Quasi Markov State Model

## Installation:

Requires Lapack, BLAS and fortran libraires:

g++ -O2 -L/usr/local/gfortran/lib -L/path_to_lapack/lib -lgfortran -lm -lblas -llapack -Wno-write-strings -o dyson-kernel dyson-kernel.cpp

## Reference:

Siqin Cao1, Andrés Montoya-Castillo, Wei Wang, Thomas E. Markland, and Xuhui Huang, On the advantages of exploiting memory in Markov state models for biomolecular dynamics, J. Chem. Phys. 153, 014105 (2020); https://doi.org/10.1063/5.0010787
