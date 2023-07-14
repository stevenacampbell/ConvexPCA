

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