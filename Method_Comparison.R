library(Rcpp)
library(microbenchmark)

#Reset Environment
rm(list=ls())

# Import source files
sourceCpp("SupportingCode/HelperFunctions_GPCA.cpp")
source("SupportingCode/HelperFunctions_ProjectedLogPCA.R")
source("SupportingCode/GPCA_Functions.R")
source("SupportingCode/GPCA_Optimizer.R")
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

### APPROXIMATE GPCA APPROACH (Campbell and Wong)

# Illustration

sub_nd<- 9 # Diadic partition for GPCA
aQFS<-approx_QF(pvals,QFS,sub_nd) # Approximate high precision quantile functions

X<-t(aQFS) # Assign data (using right endpoint convention)
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

# Specify the number of convex pcs to find
npcs<-2

# Run CPCA function
gpcs<-GPCA(npcs,X,x0,A_main,A_lower,b,h,nthreads)
gpc1<-gpcs[,1]
gpc2<-gpcs[,2]

# Plot GPCs
sub_pvals<-pvals[2:ncol(QFS)]
p<-log2(length(sub_pvals)/(2^sub_nd))
sub_pvals<-sub_pvals[seq(2^p,ncol(QFS)-1,2^p)]
plot(sub_pvals,gpc1,type="l", col="red", 
     ylim=c(min(gpcs),max(gpcs)), ylab = "Convex PCs", xlab = "Omega")
lines(sub_pvals,gpc2, type="l", col="blue")

## Reconstruction Error

# Store GPCs
W<-cbind(gpc1,gpc2)

# Project onto GPCs
Lambdas<-loadings_GPCs_fulldata(X,x0,W,A,b) # Get Loadings on GPCs
Fit<-projected_data_GPCs(X,x0,W,A,b) # Convert to function representation on grid
Fit_expanded<-expand_matrix(t(Fit),nd,sub_nd) # extend to fine grid

# Plot quantile function on fine grid and coarse approximation
plot(pvals[2:ncol(QFS)],QFS[1,2:ncol(QFS)],col="black",type="l",ylab="",xlab="")
lines(pvals[2:ncol(QFS)],Fit_expanded[1,],col="red")

# Compute Reconstruction Error
Target<-QFS[,2:ncol(QFS)] # Target (high precision) Quantile Functions
dx<-max(diff(pvals)) # Grid spacing
reconstruction_error(t(Fit_expanded),t(Target),dx) # Error

# Testing Function
runtime_gpca<-function(i){
  sub_nd<-i # Diadic partition for GPCA
  aQFS<-approx_QF(pvals,QFS,sub_nd) # Approximate high precision quantile functions
  
  X<-t(aQFS) # Assign data (using right endpoint convention)
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
  
  # Specify the number of convex pcs to find
  npcs<-2
  
  # Run CPCA function
  gpcs<-GPCA(npcs,X,x0,A_main,A_lower,b,h,nthreads)
  gpcs_<-t(expand_matrix(t(gpcs),nd,sub_nd)) # convert to natural scale for comparison
  
  return(gpcs_)
}

# Test Computation Time
time_res <- microbenchmark(
  runtime_gpca(4),
  times = 100
)
time_res

# Store two outputs for later comparison with the projected approach
gpca_<-runtime_gpca(7)
gpca_2<-runtime_gpca(9)

### PROJECTED PCA APPROACH (Beraha & Pegoraro)

# Illustration

# Grid for Quantile Functions
eval_grid<-pvals
dx <- 1/2^nd

# Create quadratic spline basis with df_-1 knot points on eval_grid
df_<-2^5
spline_basis <- bs(eval_grid,df=df_,degree=2, intercept=TRUE)

# Illustrate spline basis
color_palette <- colorRampPalette(c("blue", "red"))(df_)
plot(eval_grid, spline_basis[,1], ylim=c(0,max(spline_basis)), 
     type='l', lwd=2, col=color_palette[1], xlab="", ylab="",
     main="Quadratic B-Spline Basis Functions on [0,1]")
for(i in 2:df_){
  lines(eval_grid, spline_basis[,i], lwd=2, col=color_palette[i])
}

# Store data
X<-QFS

# Get spline coefficients
A<-fit_monotone_splines(X,spline_basis,df_)

# Obtain frechet mean
a0<-frechet_mean_coeffs(A)

# Obtain "unprojected" principal components
eigen<-solve_eigenvectors(A,a0,spline_basis,df_,dx)

# Project data onto principal components

n_ppcs<-2 # Number of PCs to extract
W<-Real(as.matrix(eigen[,1:n_ppcs])) # PCs

id<-100 # Data id for projection
a<-A[id,] # Data point

# Perform projection
E<-norm_matrix_spline(spline_basis,df_,dx)
G<-spline_constraint_matrix(df_)
v<-constr_proj_components(a,a0,W,G,E)

# Plot original data, spline fit, and projected PC fit
plot(eval_grid, X[id,],type="l",xlab="",ylab="Values",lwd=2)
lines(eval_grid, spline_basis %*% A[id,],col="blue",lwd=2)
lines(eval_grid, spline_basis %*% v,col="red",lwd="2")

# Plot principal components
pc_id <- 2 # principal component number
pc_fn_repr <- spline_basis %*% Real(eigen[,pc_id])
pc_fn_repr <- pc_fn_repr/sqrt(sum(pc_fn_repr^2)) # normalized function representation
plot(eval_grid,pc_fn_repr, type="l",col="blue",xlab="",ylab="Values")

# Reconstruction Error
Fit<-projected_data_components(A,a0,W,G,E) # Project data onto components
Fit_NS <- convert_natural_scale(spline_basis,Fit) # Fit on natural scale

# Compute Reconstruction Error
Target<-QFS[,2:ncol(QFS)] # Target (high precision) Quantile Functions
Estimate<-Fit_NS[2:ncol(QFS),] # Estimate from fits
dx<-max(diff(pvals)) # grid spacing
reconstruction_error(Estimate,t(Target),dx) # Error

# Testing Function
runtime_projected<-function(i){
  # Grid for Quantile Functions
  eval_grid<-pvals
  dx <- 1/2^nd
  
  # Create quadratic spline basis with df_-1 knot points on eval_grid
  df_<-2^i
  spline_basis <- bs(eval_grid,df=df_,degree=2, intercept=TRUE)
  
  # Store data
  X<-QFS
  
  # Get spline coefficients
  A<-fit_monotone_splines(X,spline_basis,df_)
  
  # Obtain frechet mean
  a0<-frechet_mean_coeffs(A)
  
  # Obtain "unprojected" principal components
  eigen<-solve_eigenvectors(A,a0,spline_basis,df_,dx)
  
  # Project data onto principal components
  
  n_ppcs<-2 # Number of PCs to extract
  W<-as.matrix(Real(eigen[,1:n_ppcs])) # PCs
  
  res<-spline_basis %*% W # PCs in natural scale
  return(res)
}

# Computational Time
time_res <- microbenchmark(
  runtime_projected(6),
  times = 100
)
time_res

# Store an output for later comparison
res<-runtime_projected(4)

# Comparison Plot of PCS
id<-1
plot(pvals[2:length(pvals)],gpca_[,id]/sqrt(sum(gpca_[,id]^2)),type="l",col="blue",ylab="", xlab="",main="",ylim=c(-0.055,0.045))
lines(pvals,-res[,id]/sqrt(sum(res[,id]^2)), col="red", type="l",lwd=2)
lines(pvals[2:length(pvals)],-gpca_2[,id]/sqrt(sum(gpca_2[,id]^2)),lwd=2)
lines(pvals[2:length(pvals)],gpca_[,id]/sqrt(sum(gpca_[,id]^2)),col="blue",lwd=2)
legend("bottomleft",legend=c("Fine Grid GPC", "Coarse Grid GPC", "Low Dim Projected PC"),
       col=c("black", "blue", "red"),lty=c(1,1,1),lwd=c(2,2,2))

id<-2
plot(pvals[2:length(pvals)],gpca_[,id]/sqrt(sum(gpca_[,id]^2)),type="l",col="blue",ylab="", xlab="",main="",ylim=c(0,0.04))
lines(pvals,-res[,id]/sqrt(sum(res[,id]^2)), col="red", type="l",lwd="2")
lines(pvals[2:length(pvals)],-gpca_2[,id]/sqrt(sum(gpca_2[,id]^2)),lwd="2")
lines(pvals[2:length(pvals)],gpca_[,id]/sqrt(sum(gpca_[,id]^2)),col="blue",lwd="2")
legend("bottomleft",legend=c("Fine Grid GPC", "Coarse Grid GPC", "Low Dim Projected PC"),
       col=c("black", "blue", "red"),lty=c(1,1,1),lwd=c(2,2,2))

