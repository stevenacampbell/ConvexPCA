library(nloptr)
library(Matrix)
library(profvis)
library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
library(pracma)
library(parallel)
library(readr)
library(dplyr)
library(rbenchmark)
library(matlib)

# Import rcpp helper functions
sourceCpp("HelperFunctions.cpp")

# Create functions of theta for optimization that incorporate 
# the defining input parameters (e.g. the data matrix X).

eval_CPCA_obj_closure <- function(theta,B,x0,A,b,X) {
  function(theta) eval_CPCA_obj(theta,B,x0,A,b,X)
}

eval_grad_CPCA_obj_parallel_closure <- function(theta,h,B,x0,A,b,X,nthreads) {
  function(theta) eval_grad_CPCA_obj_parallel(theta,h,B,x0,A,b,X,nthreads)
}

CPCA<-function(npcs,X,x0,A,b,h,nthreads){
  # Computes the Convex Principal Components
  
  # npcs number of principal components requested
  # X data matrix
  # x0 reference element
  # A constraint matrix
  # b constraint vector
  # h central difference step size
  # nthreads number of threads requested
  
  # Extract Dimension
  d<-nrow(X)
  
  # Allocate Matrix
  PCs<-matrix(0,ncol=npcs,nrow=d)
  
  # Optimization specifications
  opts <- list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8)
  
  for(i in seq(1,npcs,1)){
    # Get principle component estimate, projected data, and Basis Matrix
    if(i>1){
      X1<-project_orthogonal_complement(X-x0,as.matrix(PCs[,1:(i-1)]))
      B<-pracma::nullspace(t(PCs[,1:(i-1)]))
      pca <- prcomp(t(X1),center=FALSE,scale=FALSE)
      p<-pca$rotation[,1]
      p0<-ls_fit(B,p)
    }else{
      B<-eye(d)
      pca <- prcomp(t(X-x0),center=FALSE,scale=FALSE)
      p<-pca$rotation[,1]
      p0<-p
    }
    
    # Get implied theta
    theta0 <- as.vector(spherical_coords(p0))
    
    # Get functions for optimization routine
    eval_f<-eval_CPCA_obj_closure(theta,B,x0,A,b,X)
    eval_grad_f <- eval_grad_CPCA_obj_parallel_closure(theta,h,B,x0,A,b,X,nthreads)
    
    status <- paste0("Optimization Running for PC: ", i, " ...")
    print(status)
    
    # Run optimization
    res <- nloptr(
      x0          = theta0,
      eval_f      = eval_f,
      eval_grad_f = eval_grad_f,
      opts        = opts
    )
    theta_sol<-res$solution
    status <- paste0("Optimization for PC: ", i, " Complete.")
    print(status)
    
    # Store solution
    PCs[,i]<-B%*%undo_spherical_coords(theta_sol)
  }
  
  return(PCs)
}