# nri
N-way random indexing

Implementation of the N-way random indexing algorithm described in F. Sandin, B. Emruli, M. Sahlgren (2016), Random indexing of multidimensional data, Knowledge and Information Systems; http://dx.doi.org/10.1007/s10115-016-1012-2 (open access).

This code requires mtrand.h implementing the Mersenne Twister algorithm for generating pseudorandom numbers.
For further information and download see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html

The Matlab_MEX.zip archive includes the first/original implementation of the N-way random indexing algorithm,
which was used for the numerical experiments included in the manuscript and article presenting the algorithm.

The C++ implementation in randomindex.cc/h is a more recent and simplified implementation,
which includes an efficient approximation of inner products described in the paper.

![DOI: 10.5281/zenodo.53766](https://zenodo.org/badge/22352/fresan/nri.svg)
