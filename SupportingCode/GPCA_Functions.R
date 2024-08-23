
DyadicPartition <- function(a,b,n){
  # Generates the nth Dyadic partition of the interval [a,b]
  
  D <- a+ ((b-a) * ((1/2)^n)) * seq(0,2^n,1)
  
  return(D)
}

constraint_matrix_GPCA <- function(n){
  # Generates the constraint matrix for the (discrete) GPCA problem
  # defined on the nth dyadic partition of a domain [a,b]
  
  A <- matrix(0, nrow=2^n+1,ncol=2^n)
  for(i in seq(1,2^n+1,1)){
    if(i > 1){
      if(i<2^n+1){
        A[i,i] <- 1
        A[i,i-1]<--1
      }else{
        A[i,i-1] <- -1
      }
    }else{
      A[i,i] <- 1
    }
  }
  return(A)
}

constraint_vector_GPCA <-function(a,b,n){
  # Generates the constraint vector for the (discrete) GPCA problem
  # defined on the nth dyadic partition of a domain [a,b]
  
  v<-c(a,rep(0,2^n-1),-b)
  return(v)
}

GetCDFs<-function(Y,grid){
  # Input matrix Y with data distributions 
  # (samples as columns, distributions indexed by rows)
  # Input grid to evaluate CDF
  
  n_dist <- nrow(Y)
  dim <- length(grid)
  CDFs <- matrix(0,nrow=n_dist,ncol=dim)
  for(i in seq(1,n_dist,1)){
    ecdf_x <- ecdf(Y[i,])
    for(j in seq(1,dim,1)){
      CDFs[i,j]<-ecdf_x(grid[j])
    }
  }
  
  return(CDFs)
}

GetQuantiles<-function(Y,pvals){
  # Input matrix Y with data distributions 
  # (samples as columns, distributions indexed by rows)
  # Input pvals to evaluate Quantile Function
  
  n_dist <- nrow(Y)
  dim <- length(pvals)
  QFs <- matrix(0,nrow=n_dist,ncol=dim)
  for(i in seq(1,n_dist,1)){
    for(j in seq(1,dim,1)){
      QFs[i,j]<-quantile(Y[i,], probs = pvals[j])
    }
  }
  
  return(QFs)
}

loadings_GPCs<-function(z,x0,W,A,b){
  # Get loadings of a data point on the GPCs
  C_ <- 2 * eye(ncol(W))
  d_<- as.vector(- 2 * t(z-x0) %*% W)
  A_ <- - A %*% W
  b_ <- as.vector(- b + A %*% x0)
  
  sol<-quadprog(C_, d_, A = A_ , b = b_, Aeq = NULL, beq = NULL,
                lb = NULL, ub = NULL)
  return(sol$xmin)
}

loadings_GPCs_fulldata<-function(X,x0,W,A,b){
  # Get loadings of all of the data on the GPCs
  dim<-ncol(W)
  Lambdas<-matrix(0,ncol=ncol(X),nrow=dim)
  for(i in seq(1,ncol(X),1)){
    z<-X[,i]
    lambda<-loadings_GPCs(z,x0,W,A,b)
    Lambdas[,i]<-as.vector(lambda)
  }
  return(Lambdas)
}

projected_data_GPCs<-function(X,x0,W,A,b){
  # Convert the loadings to function representations of the
  # projected data.
  Lambdas<-loadings_GPCs_fulldata(X,x0,W,A,b)
  Fit<- x0 + W %*% Lambdas
  return(Fit)
}

expand_matrix<-function(M,nd,sub_nd){
  # Extend function on a Dyadic grid of size sub_nd to size nd.
  n<-nd-sub_nd
  new_matrix <- t(apply(M, 1, function(row) rep(row, each=2^n)))
  return(new_matrix)
}
