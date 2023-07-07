library(Matrix)

ilr<-function(p){
  # Isometric Log Ratio Transform
  
  dim <- length(p)
  lp <- log(p)
  cs_lp<-cumsum(lp)
  ids<-seq(1,dim,1)
  factor<-sqrt(ids[1:(dim-1)]/ids[2:dim])
  ids_<-ids[1:(dim-1)]
  x<- factor*((1/ids_)*cs_lp[1:(dim-1)]-lp[2:dim])

  return(x)
}

ilr_inv<-function(x){
  # Inverse Isometric Log Ratio Transform
  
  dim <- length(x)+1
  ids <- seq(1,dim,1)
  factor <- ids[1:(dim-1)]/ids[2:dim]
  temp <- x/sqrt(factor)
  temp[2:(dim-1)] <- temp[2:(dim-1)]-factor[1:(dim-2)]*temp[1:(dim-2)]
  cs <- cumsum(-temp)
  temp <- sum(exp(cs))
  p <- rep(0,dim)
  p[1] <- 1/(1+temp)
  p[2:dim] <- exp(cs+log(p[1]))
  
  return(p)
}

generator_matrix_ordered_simplex<-function(n){
  # Outputs the generator matrix of the n-dimensional ordered simplex
  
  dim<-n+1
  G<-matrix(0,nrow=dim-1,ncol=dim-1)
  x<-matrix(0,nrow=dim,ncol=1)
  for(i in seq(1,dim-1,1)){
    x[1:i]<-2/(dim+i)
    x[(i+1):dim]<-1/(dim+i)
    G[,i]<-ilr(x)
  }
  return(G)
}

constraint_matrix_ordered_simplex<-function(n){
  # Outputs the CPCA constraint matrix for the n-dimensional ordered simplex
  
  Ad<-matrix(0,ncol=1,nrow=n)
  Ad_<-matrix(0,ncol=1,nrow=n)
  for(k in seq(1,n,1)){
    Ad[k]<-sqrt((k+1)/k)
    if(k<n){
      Ad_[k+1]<--sqrt(k/(k+1))
    }
  }
  Ad<-Ad/log(2)
  Ad_<-Ad_[2:length(Ad_)]/log(2)
  A<-bandSparse(n, n, (-1):0, list(Ad_, Ad))
  
  return(as.matrix(A))
}


#a<-c(0.2,0.7,0.1)
#ilr(a)
#ilr_inv(ilr(a))

#err<-rep(0,1000)
#for(i in seq(1,1000,1)){
#  x<-runif(10)
#  y<-ilr_inv(x)
#  err[i]<-sum(abs(x-ilr(y)))+sum(abs((ilr_inv(ilr(y))-y)))
#}
#plot(err)
