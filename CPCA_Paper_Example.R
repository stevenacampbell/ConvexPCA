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
library("gg3D")
library(tidyverse)
library(rbenchmark)
library(matlib)

#Reset Environment
rm(list=ls())

#Problem Setup

# Let d be the dimension of the original problem in R^d.
# After the first j principle components found we find a matrix B of orthogonal basis vectors with dimension d x (d-j)
# We introduce a new vector theta for this problem of dimension (d-j-1) 
# This corresponds to a vector p(theta) of dimension (d-j) which is multiplied with B to get the input for our problem
# The new vector Bp(theta) is a linear combination of the orthogonal basis vectors and will be orthogonal to the earlier PCs
# The problem then will return a new minimum in terms of theta* from which we can recover the new PC Bp(theta*)
# We repeat this process where of course the first principle component is found by setting B=Id
# We can find candidate solutions near the local minimum by using the usual PCA solution p_ for the projected data (see notes).
# In order to get our initial theta0 from this we must solve p_=Bp(theta) for theta

#Least Squares Fit to solve for solution to system Bp=p_
LSFit<-function(B,b){
  L<-t(B)%*%B
  d<-t(B)%*%b
  p<-solve(L,d)
  return(p)
}

#Create Matrix B
GetB<-function(PCs,d,npcs){
  m<-d-npcs
  B<-eye(d)
  B<-B[,1:m]
  B1<-B
  for(i in seq(1,ncol(B),1)){
    if(npcs>1){
      for(j in seq(1,npcs,1)){
        B1[,i]<-B1[,i]-c((B[,i]%*%PCs[,j])/(PCs[,j]%*%PCs[,j]))*PCs[,j]
      }
    }else{
      B1[,i]<-B1[,i]-c((B[,i]%*%PCs)/(PCs%*%PCs))*PCs
    }
  }
  B<-GramSchmidt(B1)
  return(B)
}

DataProjection<-function(X0,X,PCs,npcs){
  nobs<-ncol(X)
  B<-X-X0
  B1<-B
  for(i in seq(1,nobs,1)){
    if(npcs>1){
      for(j in seq(1,npcs,1)){
        B1[,i]<-B1[,i]-c((B[,i]%*%PCs[,j])/(PCs[,j]%*%PCs[,j]))*PCs[,j]
      }
    }else{
      B1[,i]<-B1[,i]-c((B[,i]%*%PCs)/(PCs%*%PCs))*PCs
    }
  }
  B<-B1
  return(B)
}

InitialP<-function(B,p){
  P<-LSFit(B,p)
  return(P)
}


################################################
################ CPCA FUNCTIONS ################
################################################

#Import CPCA Cpp Functions
setwd("~/CapitalDistributionProject")
sourceCpp("ParCPCATest.cpp")


#ILR TRANSFORM
ILR <-function(p){
  #p has observations as columns and variables as rows
  x<-matrix(0,ncol=ncol(p),nrow=nrow(p)-1)
  for(j in seq(1,ncol(p),1)){
    for(i in seq(1,nrow(p)-1,1)){
      x[i,j]<-sqrt(i/(i+1))*log((exp(mean(log(p[1:i,j]))))/p[i+1,j])
    }
  }
  return(x)
}

##INVERSE ILR TRANSFORM
ILR_INV <- function(x){
  #x has observations as columns and variables as rows
  #returns vector in the simplex
  II<-as.matrix(rep(0,nrow(x)))
  for(i in seq(1,nrow(x),1)){
    II[i]<-sqrt((i+1)/i)
  }
  b<-x
  for(i in seq(1,ncol(b),1)){
    b[,i]<-exp(b[,i]*II)
  }
  c<-b
  for(i in seq(1,nrow(b)-1,1)){
    b[i+1,]<-(c[i+1,]/(c[i,])^(i/(i+1)))
  }
  phat<-matrix(1,ncol=ncol(b),nrow=nrow(b)+1)
  for(i in seq(1,nrow(b),1)){
    phat[i+1,]<-phat[i,]/b[i,]
  }
  for(j in seq(1,ncol(phat),1)){
    phat[,j]<-phat[,j]/sum(phat[,j])
  }
  return(phat)
}

#GENERATOR MATRIX:
Generator<-function(dim){
  G<-matrix(0,nrow=dim-1,ncol=dim-1)
  x<-matrix(0,nrow=dim,ncol=1)
  for(i in seq(1,dim-1,1)){
    x[1:i]<-2/(dim+i)
    x[(i+1):dim]<-1/(dim+i)
    xILR<-ILR(cbind(x,x))
    G[,i]<-xILR[,1]
  }
  return(G)
}

#Matrix diagonals
AVECTORS<-function(dim){
  Ad<-matrix(0,ncol=1,nrow=dim-1)
  Ad_<-matrix(0,ncol=1,nrow=dim-1)
  for(k in seq(1,dim-1,1)){
    Ad[k]<-sqrt((k+1)/k)
    if(k<dim-1){
      Ad_[k+1]<--sqrt(k/(k+1))
    }
  }
  Ad<-Ad/log(2)
  Ad_<-Ad_[2:length(Ad_)]/log(2)
  A<-bandSparse(dim-1, dim-1, (-1):0, list(Ad_, Ad))
  return(list(Ad,Ad_,A))
}

#NUMERATOR FOR INTERSECTION POINT
Numerator<-function(X0,A){
  NUM<-as.vector(A %*% X0)
  return(NUM)
}

#Recover Spherical Coords from Vector
SPHERICAL<-function(P){
  if(length(P)==2){
    theta=acos(P[1])
  }else{
    P2<-P*P
    P2r<-rev(P2)
    CS<-cumsum(P2r)
    CSr<-rev(CS)
    d<-length(P)
    CSr_<-CSr[1:(d-2)]
    P_<-P[1:(d-2)]
    theta<-acos(P_/sqrt(CSr_))
    if(P[d]<0){
      v<-2*pi-acos(P[d-1]/sqrt(CSr[d-1]))
    }else{
      v<-acos(P[d-1]/sqrt(CSr[d-1]))
    }
    theta<-as.vector(c(theta,v))
  }
  return(theta)
}

eval_f_par<-function(theta,B,len){
  FunctionEval(B,Ad,Ad_,X0,theta,NUM,X,nobs,inobs,d,len)
}

eval_grad_f_par<-function(theta,B,len){
  GradientEval(B,Ad,Ad_,X0,theta,NUM,X,nobs,inobs,d,len,h)
}

CPCA<-function(X,npcs,len,d){
  PCs<-matrix(0,ncol=npcs,nrow=d)
  StdPCs<-matrix(0,ncol=npcs,nrow=d)
  
  X0 <-as.vector(rowMeans(X))
  
  opts <- list("algorithm"="NLOPT_LD_LBFGS",
               "xtol_rel"=1.0e-8)
  for(i in seq(1,npcs,1)){
    
    #Get principle component estimate, projected data, and Basis Matrix alongside P0
    if(i>1){
      X1<-DataProjection(X0,X,PCs[,1:(i-1)],i-1)
      B<-GetB(PCs[,1:(i-1)],d,i-1)
      pca <- prcomp(t(X1),center=FALSE,scale=FALSE)
      P<-pca$rotation[,1]
      StdPCs[,i]<-pca$rotation[,1]
      P0<-InitialP(B,P)
    }else{
      B<-eye(d)
      pca <- prcomp(t(X-X0),center=FALSE,scale=FALSE)
      P0<-pca$rotation[,1]
      StdPCs[,i]<-pca$rotation[,1]
    }
    
    #GET IMPLIED THETA
    theta0 <- as.vector(SPHERICAL(P0))
    
    l0 = len-(i-1)
    
    status <- paste0("Optimization Running for PC: ", i, " ...")
    print(status)
    res <- nloptr(
      x0          = theta0,
      eval_f      = eval_f_par,
      eval_grad_f = eval_grad_f_par,
      opts        = opts,
      B=B,
      len=l0)
    theta_sol<-res$solution
    status <- paste0("Optimization for PC: ", i, " Complete.")
    print(status)
    
    PCs[,i]<-B%*%GetP_v(theta_sol,len-(i-1))
  }
  
  return(list(PCs,StdPCs))
}



#########################################################
####################### CPCA TEST #######################
#########################################################
#GLOBAL VARIABLES
setThreadOptions(numThreads = 4) #Number of Threads
#defaultNumThreads()
dim<-3
d<-dim-1
len<-d-1
zero<-matrix(0,nrow=dim-1,ncol=1)
ZERO<-matrix(0,nrow=dim-1,ncol=dim-1)
ids<-matrix(-1,nrow=2,ncol=1)
h<-0.000001

set.seed(100)

#GET MATRIX
Alist<-AVECTORS(dim)
Ad<-Alist[[1]]
Ad_<-Alist[[2]]
A<-Alist[[3]]

#TEST SAMPLE
num<-1000
P<-matrix(0,ncol=num,nrow=dim)
for(i in seq(1,num,1)){
  #Generate random distribution
  #generate random vector of size one less than dim
  rand<-runif(dim-1)
  #sort these from least to greatest
  rand<-rand[order(rand)]
  #append 0 and 1 to the end of this vector and compute the differences
  #this is an independent uniform sample from the unit simplex
  rand<-c(0,rand,1)
  r<-diff(rand)
  P[,i]<--sort(-r)
}
X<-ILR(P)
nobs<-ncol(X)
inobs<-1/nobs
X0<-as.vector(rowMeans(X))

G<-Generator(dim)

X1<-X
mask<-(X1[1,]>0.15)
X1<-X1[,mask]
mask<-(X1[1,]<1.0)
X1<-X1[,mask]
X0<-as.vector(rowMeans(X1))
plot(0:5, 0:5, xlab="x",ylab="y",type = "n",main="Projections")# setting up coord. system
abline(0,G[2,1]/G[1,1], lwd=3, col="black")
abline(v=0,lwd=3, col="black")
points(X1[1,],X1[2,], col="blue", lwd=10,pch=17)
points(X0[1],X0[2],col="red",lwd=10,pch=15)
abline()

#GET NUMERATOR VALUE FOR INTERSECTION
NUM<-as.vector(Numerator(X0,A))

#Get principle component estimate
pca <- prcomp(t(X-X0),center=FALSE,scale=FALSE)
summary(pca)
P_<-pca$rotation[,1]

B<-eye(d)

#GET IMPLIED THETA
theta0 <- as.vector(SPHERICAL(P_))

#SOLVE CPCA PROBLEM
npcs<-1

sol<-CPCA(X1,npcs,len,d)
CPCs<-sol[[1]]
StdPCs<-sol[[2]]

plot(seq(0,1.65,0.165), seq(0,5,0.5), xaxt = "n",yaxt = "n", xlab="", ylab="", type = "n")
clip(0,1.695, 0, 5.1)
polygon(c(seq(0,1.8,0.18),seq(1.8,0,-0.18)),c(rep(6,11),(G[2,1]/G[1,1])*seq(1.8,0,-0.18)), col = "grey95")
abline(0,G[2,1]/G[1,1],lwd=3)
abline(v=0,lwd=5)
points(X1[1,],X1[2,],pch=0, lwd=1
       ,col="grey50")
P<-CPCs[,1]
m<-(X0[2]+P[2]-X0[2])/(X0[1]+P[1]-X0[1])
b<-X0[2]-m*X0[1]
z<-seq(0,1.65,0.00165)
lines(z[190:length(z)],m*z[190:length(z)]+b,col="blue",lwd=3)
points(X0[1],X0[2],col="red",pch=15,lwd=10)
rowMeans(X)
abline(0,G[2,1]/G[1,1],lwd=3)
abline(v=0,lwd=5)
pca <- prcomp(t(X1-X0),center=FALSE,scale=FALSE)
P_<-pca$rotation[,1]
m<-(X0[2]+P_[2]-X0[2])/(X0[1]+P_[1]-X0[1])
b<-X0[2]-m*X0[1]
for(k in seq(-3,3,0.5)){
  abline(b+k*0.01,m, col="purple",lwd=2, lty=2)
}
points(X0[1],X0[2],col="red",pch=15,lwd=5)
legend( 1.1255,1.551,#x="bottomright", 
        legend=c("Convex PC","Euclidean PC","Convex Set Bdy","Data Points", "Reference Point"), 
        col=c("blue","purple","black","grey50","red"), lty=c(1,2,1,NA,NA), 
        pch=c(NA,NA,NA,0,15), lwd=c(2.9,2.9,2.9,1,10), merge=FALSE )
