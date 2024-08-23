library(Rcpp)
library(pracma)
library(microbenchmark)

#Reset Environment
rm(list=ls())

# Import source files
sourceCpp("SupportingCode/HelperFunctions_GPCA.cpp")
source("SupportingCode/HelperFunctions_ProjectedLogPCA.R")
source("SupportingCode/GPCA_Functions.R")
source("SupportingCode/GPCA_Optimizer.R")
source("SupportingCode/AitchisonGeometryOperations.R")
source("SupportingCode/AtlasModelFunctions.R")
source("SupportingCode/Misc.R")

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

# Obtain Empirical Quantile Functions (High Precision)
nd<-12 # nth Dyadic partition for Omega grid
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

### Test Function Evaluation Times for Approximate GPCA Approach (Campbell and Wong)

sub_nd<- 7 # Diadic partition for GPCA (Vary this dimension for computational times)
aQFS<-approx_QF(pvals,QFS,sub_nd) # Approximate high precision quantile functions

X<-t(aQFS) # Assign data
x0<-rowMeans(X) # compute mean distribution

h<-1e-10 # central difference step size

# Get constraints
A<-constraint_matrix_GPCA(sub_nd)
b<-constraint_vector_GPCA(omega0,omega1,sub_nd)

# Extract Diagonals
A_main<-diag(A)
A_lower<-lower_off_diagonal(A)

# Specify the number of threads to be used
nthreads<-1

# Generate random points for function evaluation
p<-runif(length(x0)) # random point for function in natural scale
theta<-spherical_coords(p) # spherical coordinates
B<-eye(length(p)) # basis matrix

# Objective Function Evaluation
time_res_obj <- microbenchmark(
  eval_GPCA_obj(p,x0,A_main,A_lower,b,X),
  times = 1000
)
time_res_obj

# Objective Function Evaluation in Spherical Coordinates
time_res_obj_sph <- microbenchmark(
  eval_GPCA_obj_sph(theta,B,x0,A_main,A_lower,b,X),
  times = 1000
)
time_res_obj_sph

# Gradient Evaluation
time_res_grad <- microbenchmark(
  eval_grad_GPCA_obj_parallel(p,h,x0,A_main,A_lower,b,X,nthreads),
  times = 1000
)
time_res_grad

# Gradient Evaluation in Spherical Coordinates
time_res_grad_sph <- microbenchmark(
  eval_grad_GPCA_obj_sph_parallel(theta,B,h,x0,A_main,A_lower,b,X,nthreads),
  times = 1000
)
time_res_grad_sph

# Test repeated inner product evaluations

# Function for redundant inner products (Note: n is the # of securities)
test_eval<-function(p,q){
  g<-0
  for(i in seq(1,n,1)){
    g<-p %*% q
  }
  return(g)
}

dimen <- 2^sub_nd # Data dimension
p<-runif(dimen) # Random point
q<-runif(dimen) # Random point

time_res_inner_prod <- microbenchmark(
  test_eval(p,q),
  times = 1000
)
time_res_inner_prod
