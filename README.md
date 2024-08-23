# Code supplement for ``Convex PCA with applications to Wasserstein geodesic PCA and ranked data"

This repository contains sample code for an implementation of Convex PCA as seen in [Campbell & Wong (2022)]. Sample applications to capital distribution curves and return distributions by rank arising from the (generalized) Atlas model of equity markets (see e.g. [Banner, Fernholz, and Karatzas (2005)]) are provided for illustration purposes. Code reproducing the simple 2d example from Figure 1 in [Campbell & Wong (2022)] is also included. All code required to reproduce the numerical experiments of Appendix B.4 in [Campbell & Wong (2022)] is provided. These experiments compare the present methodology with the approaches of [Cazelles et al. (2018)] and [Pegoraro and Beraha (2022)] for GPCA and the dimensionality reduction of distributional data. A more detailed breakdown of the files and their contents is included here:

- SupportingCode Folder: Folder containing all supporting function definitions in R and C++.
- ToyExample_2D.R: Implementation of the 2d example from Figure 1 in [Campbell & Wong (2022)].
- ToyExample_CPCA_CapitalDistributionCurves_AtlasModel.R: Applications to capital distribution curves.
- ToyExample_GPCA_RankReturnDistributions_AtlasModel.R: Applications to return distributions by rank.
- Method_Comparison.R: Numerical comparison with the approach of [Pegoraro and Beraha (2022)].
- FunctionEvaluation_ComputationTime.R: Test of the computational cost of function and gradient evaluations for GPCA.
- DataPrep_Cazellesetal_Comparison.R: Preparation of synthetic data for a comparison with the approach of [Cazelles et al. (2018)]. This data is used as an input to their publicly available MATLAB code (see https://github.com/ecazelles/2017-GPCA-vs-LogPCA-Wasserstein/tree/master). 
- CazellesTest.m: MATLAB File that feeds the histogram and bin data into the functions of [Cazelles et al. (2018)]. PLEASE NOTE: Their code must be downloaded separately and stored in the same directory in order for this file to run.


NOTE: The CPCA functionality included here allows for arbitrary constraint matrix-vector pairs (A,b) and is implemented in C++ using RCPP. The optimization routine allows for parallelization. The number of threads to be recruited for this task is specified by the user using the parameter "nthreads" in the example files. Additional efficiency gains may be obtained if a particular structure for A and b is assumed (e.g. if A is a tri-diagonal matrix) and should be considered on a case-by-case basis. The provided implementation of GPCA incorporates these considerations for the bi-diagonal constraint matrix that arises in that context.
