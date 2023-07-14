library(Rcpp)
library(RcppParallel)

#Reset Environment
rm(list=ls())

# Import source files
sourceCpp("HelperFunctions.cpp")
source("CPCAFunctions.R")
source("AitchisonGeometryOperations.R")
source("AtlasModelFunctions.R")
source("GPCAFunctions.R")
source("Misc.R")

# Set seed
set.seed(100)

# Initialize Parameters for Atlas Model
n<-101 # number of securities
N<-10000 # number of time steps
T<-100 # terminal time
gamma<-0.1 # Atlas baseline drift parameter
gs<-((-seq(n,1,-1)+(n+1)/2)/((n+1)/2)) # Atlas rank based drifts
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

# Compute the returns and returns by rank
r<-Returns(M)
rbr<-ReturnsbyRank(M,r)


### Data Illustration ###

# Define range of data
omega0 <- min(rbr)
omega1 <- max(rbr)

# Define the breaks for the return histograms
breaks <- seq(omega0, omega1, length.out = 50)

# Calculate rank based return histograms
hist_1 <- hist(rbr[1,], breaks = breaks, plot = FALSE)
hist_2 <- hist(rbr[26,], breaks = breaks, plot = FALSE)
hist_3 <- hist(rbr[51,], breaks = breaks, plot = FALSE)
hist_4 <- hist(rbr[76,], breaks = breaks, plot = FALSE)
hist_5 <- hist(rbr[101,], breaks = breaks, plot = FALSE)

# Plot rank return histograms
color_palette <- colorRampPalette(c("blue","red"))(5)
c1<-transparent_col(color_palette[1],0.8) # Adjust color transparency
c2<-transparent_col(color_palette[2],0.8)
c3<-transparent_col(color_palette[3],0.8)
c4<-transparent_col(color_palette[4],0.8)
c5<-transparent_col(color_palette[5],0.8)
plot(hist_1, col = c1, main = "Rank Return Distributions", 
     xlab = "Returns", ylab = "Frequency")
plot(hist_2, col = c2, add = TRUE)
plot(hist_3, col = c3, add = TRUE)
plot(hist_4, col = c4, add = TRUE)
plot(hist_5, col = c5, add = TRUE)
legend("topright", legend = c(1,26,51,76,101), fill = c(c1,c2,c3,c4,c5))

# Obtain Empirical CDF
grid<-seq(omega0,omega1,0.01) # Define grid for CDF values
CDFS<- GetCDFs(rbr,grid) # Obtain CDF values on the grid

# Obtain Empirical Quantile Functions
nd<-8 # nth Dyadic partition for Omega grid
omega_grid <- DyadicPartition(omega0,omega1,nd) # Obtain dyadic parition
pvals <- punif(omega_grid, min=omega0, max=omega1) # Convert to p values using uniform reference measure
QFS<- GetQuantiles(rbr,pvals) # Obtain Quantile function values

# Plot CDFs
color_palette <- colorRampPalette(c("blue","red"))(n)
plot(grid,CDFS[1,],type="l", col = color_palette[1], ylim=c(min(CDFS),max(CDFS)))
for(i in seq(1,n,1)){
  lines(grid,CDFS[i,],col=color_palette[i])
}

# Plot Quantile Functions
color_palette <- colorRampPalette(c("blue","red"))(n)
plot(pvals,QFS[1,],type="l", col = color_palette[1], ylim=c(min(QFS),max(QFS)))
for(i in seq(1,n,1)){
  lines(pvals,QFS[i,],col=color_palette[i])
}

# Plot Transport function 
# (by using the uniform distribution as reference measure)
color_palette <- colorRampPalette(c("blue","red"))(n)
plot(omega_grid,QFS[1,],type="l", col = color_palette[1], ylim=c(min(QFS),max(QFS)))
for(i in seq(1,n,1)){
  lines(omega_grid,QFS[i,],col=color_palette[i])
}

### CPCA Analysis Set Up ###

X<-t(QFS[,2:ncol(QFS)]) # Assign data (using right endpoint convention)
x0<-rowMeans(X) # compute mean distribution

h<-1e-10 # central difference step size

# Get constraints
A<-constraint_matrix_GPCA(nd)
b<-constraint_vector_GPCA(omega0,omega1,nd)

# Specify the number of threads to be used
nthreads<-8

# Specify the number of convex pcs to find
npcs<-2

# Run CPCA function

# Note: The dyadic partition allows us to directly use the (discretized) 
# transport maps in the CPCA function without re-weighting as the induced 
# distance is the same up to a constant.

cpcs<-CPCA(npcs,X,x0,A,b,h,nthreads)
cpc1<-cpcs[,1]
cpc2<-cpcs[,2]

# Plot Convex PCs
plot(omega_grid[2:length(pvals)],cpc1,type="l", col="red", 
     ylim=c(min(cpcs),max(cpcs)), ylab = "Convex PCs", xlab = "Omega")
lines(omega_grid[2:length(pvals)],cpc2, type="l", col="blue")


### QUANTILE FIGURE - 1st CPC ###

epsilons<- seq(0,1,0.01) # specify set of epsilons for perturbations

# Initialize perturbation matrix to store perturbed quantile functions
perturbed_plus<-matrix(0,nrow=length(epsilons),ncol=length(x0))
perturbed_minus<-matrix(0,nrow=length(epsilons),ncol=length(x0))

# Choose which convex principal direction to perturb along
pertdir<-1

# Compute perturbed quantile functions using the data mean as a reference point
index<-1
for(epsilon in epsilons){
  for(i in seq(0,1,1)){
    if(i==0){
      x_m<-x0-epsilon*cpcs[,pertdir]
    }else{
      x_p<-x0+epsilon*cpcs[,pertdir]
    }
  }
  perturbed_plus[index,]<-x_p
  perturbed_minus[index,]<-x_m
  index<-index+1
}

# Specify plot range for x axis
xmin<- -0.7 # minimum x axis value
xmax<- 0.7 # maximum x axis value

# Initialize a list of (upper) quantiles to plot
# (the result will plot 1-quantiles specified as well)
quantiles_list <- seq(0.5,0.99,0.01)

# Obtain color gradient
colfunc<-colorRampPalette(c("red","blue"))
cols<-colfunc(length(quantiles_list))

# Plot perturbed quantiles
vals <- perturbed_plus[,min(which(pvals>=0.5))]
plot(epsilons~vals, ylim=c(-max(epsilons),max(epsilons)),xlim=c(xmin,xmax),type="l",col=cols[1],ylab="t",xlab="Return")
for(i in seq(1,length(quantiles_list),1)){
  vals <-perturbed_plus[,min(which(pvals>=quantiles_list[i]))]
  lines(epsilons~vals, ylim=c(-max(epsilons),max(epsilons)),xlim=c(xmin,xmax),type="l",col=cols[i])
}
quantiles_list_ <- 1-quantiles_list
for(i in seq(1,length(quantiles_list),1)){
  vals <-perturbed_plus[,min(which(pvals>=quantiles_list_[i]))]
  lines(epsilons~vals, ylim=c(-max(epsilons),max(epsilons)),xlim=c(xmin,xmax),type="l",col=cols[i])
}
for(i in seq(1,length(quantiles_list),1)){
  vals <-perturbed_minus[,min(which(pvals>=quantiles_list[i]))]
  lines(-epsilons~vals, ylim=c(-max(epsilons),max(epsilons)),xlim=c(xmin,xmax),type="l",col=cols[i])
}
for(i in seq(1,length(quantiles_list),1)){
  vals <-perturbed_minus[,min(which(pvals>=quantiles_list_[i]))]
  lines(-epsilons~vals, ylim=c(-max(epsilons),max(epsilons)),xlim=c(xmin,xmax),type="l",col=cols[i])
}

