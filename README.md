# Online code supplement for ``Convex PCA with applications to Wasserstein geodesic PCA and ranked data"

This repository contains sample code for an implementation of Convex PCA as seen in [Campbell & Wong (2022)]. Sample applications to Capital Distribution Curves and Return Distributions by Rank arising in the (generalized) Atlas model of equity markets are provided for illustration purposes. Code reproducing the simple 2d example from Figure 1 in [Campbell & Wong (2022)] is also included.

Note: The CPCA function included here allows for arbitrary constraint matrix-vector pairs (A,b) and is implemented in C++ using RCPP. The optimization routine allows for parallelization. The number of threads to be recruited for this task is specified by the user using the parameter "nthreads" in the example files. Additional efficiency gains may be obtained if a particular structure for A and b is assumed (e.g. if A is a tri-diagonal matrix) and should be considered on a case-by-case basis.
