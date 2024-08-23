# Load necessary libraries
library(splines)
library(pracma)
library(quadprog)

spline_constraint_matrix<-function(n_funct){
  # Constraint matrix for monotone splines
  
  G<-eye(n_funct,n_funct)
  for(i in seq(1,n_funct,1)){
    G[i,i-1]<--1
  }
  G<-G[-1,]
  
  return(G)
}

fit_monotone_splines<-function(Y, spline_basis, n_funct){
  # Input: Matrix where Y has the functions as rows and values as columns,
  # spline basis and number of functions
  # Output: Matrix of coefficients of constrained least squares fit
  
  coeffs <- zeros(nrow(Y),n_funct)
  
  G<-spline_constraint_matrix(n_funct) # Constraint Matrix
  
  b<-rep(0,n_funct-1) # Constraint vector
  
  C<- 2*t(spline_basis) %*% spline_basis # Cost matrix
  
  for(i in seq(1,nrow(Y),1)){
    d <- as.vector(-2*t(Y[i,]) %*% spline_basis) # Cost vector
    
    # solve constrained quadratic program for spline coefficients
    sol<-quadprog(C, d, A = -G, b = b, Aeq = NULL, beq = NULL,
                  lb = NULL, ub = NULL)
    coeffs[i,]<-sol$xmin
  }

  return(coeffs)
}

norm_matrix_spline<-function(spline_basis,n_funct,dx){
  # Output norm matrix for splines from P & B (2022)
  E<-zeros(n_funct,n_funct)
  for(i in seq(1,n_funct,1)){
    for(j in seq(1,n_funct,1)){
      E[i,j]<-sum(spline_basis[,i]*spline_basis[,j]*dx)
    }
  }
  return(E)
}

projection_monotone_coefficients<-function(v,spline_basis,n_funct,dx){
  # Project onto montone coefficients
  
  E<-norm_matrix_spline(spline_basis,n_funct,dx)
  G<-spline_constraint_matrix(n_funct)
  b<-rep(0,n_funct-1)
  d<- as.vector(-2 * E %*% v)
  
  sol<-quadprog(2*E, d, A = -G, b = b, Aeq = NULL, beq = NULL,
                lb = NULL, ub = NULL)
  
  return(sol$xmin)
}

frechet_mean_coeffs<-function(A){
  # Compute Frechet mean in the coefficient space
  a0<-colMeans(A)
  return(a0)
}

solve_eigenvectors<-function(coeffs,a0,spline_basis,n_funct,dx){
  # Solve for eigenvectors that determine PCs
  A <- sweep(coeffs,2,a0,"-")
  E <- norm_matrix_spline(spline_basis,n_funct,dx)
  M <- t(A) %*% A %*% E
  sol <- eigen(M)
  v <- sol$vectors
  v <- sweep(v,2,sqrt(as.vector(diag(t(v) %*% E %*% v))),"/")
  return(v)
}

loadings_constr_proj_components<-function(a,a0,W,G,E){
  # Extract loadings on PCs of data
  C_ <- 2 * eye(ncol(W))
  d_<- as.vector(- 2 * t(a-a0) %*% E %*% W)
  A_ <- - G %*% W
  b_ <- as.vector(G %*% a0)
  
  sol<-quadprog(C_, d_, A = A_ , b = b_, Aeq = NULL, beq = NULL,
                lb = NULL, ub = NULL)
  return(sol$xmin)
}

constr_proj_components<-function(a,a0,W,G,E){
  # Extract coefficient representation of data projected onto the 
  # "projected" principal component from loadings
  lambda <- loadings_constr_proj_components(a,a0,W,G,E)
  
  return(a0 + W %*% lambda)
}


loadings_constr_proj_fulldata<-function(A,a0,W,G,E){
  # Extract loadings for the full data
  dim<-ncol(W)
  Lambdas<-matrix(0,ncol=nrow(A),nrow=dim)
  for(i in seq(1,nrow(A),1)){
    a<-A[i,]
    lambda<-loadings_constr_proj_components(a,a0,W,G,E)
    Lambdas[,i]<-as.vector(lambda)
  }
  return(Lambdas)
}

projected_data_components<-function(A,a0,W,G,E){
  # Extract coefficient representation of data projected onto the 
  # "projected" principal component for the full data
  Lambdas<-loadings_constr_proj_fulldata(A,a0,W,G,E)
  Fit<- a0 + W %*% Lambdas
  return(Fit)
}

convert_natural_scale<-function(spline_basis,Fit){
  # Convert spline coefficients to function representation
  Fit_NS <- spline_basis %*% Fit
  return(Fit_NS)
}

w2_approx_error <- function(q1,q2,dx){
  # Compute W2 error from quantile functions on a uniform grid 
  # with spacing dx.
  return(sqrt(sum((q1-q2)^2*dx)))
}

reconstruction_error <- function(Q,Q_,dx){
  # Compute reconstruction error from 2 sets of quantiles Q.
  # Note: columns index the observations.
  err <- 0
  n <- ncol(Q)
  for(i in seq(1,n,1)){
    err <- err + (1/n) * w2_approx_error(Q[,i],Q_[,i],dx)
  }
  return(err)
}

