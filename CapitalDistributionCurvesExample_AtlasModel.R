library(Rcpp)
library(RcppParallel)

#Reset Environment
rm(list=ls())

# Import source files
sourceCpp("HelperFunctions.cpp")
source("CPCAFunctions.R")
source("AitchisonGeometryOperations.R")
source("AtlasModelFunctions.R")

# Initialize Parameters for Atlas Model
n<-21 # number of securities
N<-10000 # number of time steps
T<-100 # terminal time
gamma<-0.1 # Atlas baseline drift parameter
gs<-((-seq(n,1,-1)+(n+1)/2)/((n+1)/2)) #Atlas rank based drifts
sigmas<-1+2*seq(0,n-1,1)/(n-1) # Atlas rank based volatilities
times<-seq(0,T,T/N) # time vector
B<-100 # burn in period
N_B<-10000 # number of time steps burn in period

# Burn in period
m0<-runif(n) #Arbitary initial market distribution
M<-AtlasModel(n,N_B,B,gamma,gs,sigmas,m0)

# Simulate System
m0<-M[,N_B] #Extract new initial distribution
M<-AtlasModel(n,N,T,gamma,gs,sigmas,m0)

# Obtain Ranked System
M_R<-RankedSystem(M)

# Plot ranked system
color_palette <- colorRampPalette(c("blue", "red"))(n)
plot(times,M_R[1,],type="l",ylim=c(min(M_R),max(M_R)),col=color_palette[1],
     ylab="Log Caps", xlab="Time")
for(i in seq(2,n,1)){
  lines(times,M_R[i,],col=color_palette[i])
}

# Compute and plot capital distribution curves
ranks<-seq(1,n,1)
color_palette <- colorRampPalette(c("blue", "red"))(N)
cap_dist <- CapitalDistribution(M)
incr <- 100 # number of time steps between plotted/sampled curves
plot(log(ranks),log(cap_dist[,1]),col=color_palette[1],
     ylim=c(min(log(cap_dist)),max(log(cap_dist))),type='l',
     ylab="Log Weights",xlab="Log Rank")
for(i in seq(2,N,incr)){
  lines(log(ranks),log(cap_dist[,i]),col=color_palette[i])
}

#Select the sampling frequency of the capital curves to be used for analysis
freq<-10
subset_cap_dist<-cap_dist[,seq(1,N,freq)]

X<-apply(subset_cap_dist,2,ilr) # compute ilr transformation of capital curves
x0<-rowMeans(X) # mean capital curve (in Aitchison geometry)
h<-1e-10 # central difference step size

# Get constraints
A<-constraint_matrix_ordered_simplex(n-1)
b<-rep(0,nrow(A))

# Specify the number of threads to be used
nthreads<-8

# Specify the number of convex pcs to find
npcs<-2

# Run CPCA function
cpcs<-CPCA(npcs,X,x0,A,b,h,nthreads)
cpc1<-cpcs[,1]
cpc2<-cpcs[,2]

# Perturb mean capital distribution
eps<-4
curves<-cbind(ilr_inv(x0-eps*cpc1),ilr_inv(x0),ilr_inv(x0+eps*cpc1))

plot(log(ranks),log(curves[,2]),ylim=c(min(log(curves)),max(log(curves))),
     type="l", ylab="Log Weights", xlab="Log Ranks")
lines(log(ranks),log(curves[,3]),type="l",col="blue")
lines(log(ranks),log(curves[,1]),type="l",col="red")

eps<-3
curves<-cbind(ilr_inv(x0-eps*cpc2),ilr_inv(x0),ilr_inv(x0+eps*cpc2))

plot(log(ranks),log(curves[,2]),ylim=c(min(log(curves)),max(log(curves))),
     type="l", ylab="Log Weights", xlab="Log Ranks")
lines(log(ranks),log(curves[,3]),type="l",col="blue")
lines(log(ranks),log(curves[,1]),type="l",col="red")

