# nri
N-way random indexing

Implementation of the N-way random indexing algorithm.

This code requires mtrand.h implementing the Mersenne Twister algorithm for generating pseudorandom numbers.
For further information and download see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html

The Matlab_MEX.zip archive includes the first/original implementation of the N-way random indexing algorithm,
which was used for the numerical experiments included in the manuscript and paper presenting this algorithm:
Fredrik Sandin, Blerim Emruli, Magnus Sahlgren (2016), Random indexing of multi-dimensional data.

The C++ implementation in randomindex.cc/h is a more recent and simplified implementation,
which includes an efficient approximation of inner products described in the paper.
