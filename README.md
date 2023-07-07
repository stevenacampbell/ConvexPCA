# ConvexPCA
Sample code for an implementation of Convex PCA. A toy application to Capital Distribution Curves arising in the (generalized) Atlas model of equity markets is provided for illustration purposes.

Note: The CPCA function included here allows for arbitrary constraint matrix-vector pairs (A,b) and is implemented in C++ using RCPP. Additional efficiency gains may be obtained if a particular structure for A and b is assumed (e.g. if A is a tri-diagonal matrix) and should be considered on a case-by-case basis.
