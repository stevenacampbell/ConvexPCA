#Reset Environment
rm(list=ls())

# Import source files
source("SupportingCode/GPCA_Functions.R")
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

nd<-7 # nth Dyadic partition for Omega grid
omega_grid <- DyadicPartition(omega0,omega1,nd)

# Construct data input for Cazelles
mu<-matrix(0,nrow=nrow(rbr),ncol=length(omega_grid)-1)
for(i in seq(1,nrow(rbr),1)){
  hist_data <- hist(rbr[i,],breaks=omega_grid, plot=FALSE)
  mu[i,]<-hist_data$density
}

# Illustrate Histograms
color_palette <- colorRampPalette(c("blue","red"))(n)
plot(omega_grid[1:(length(omega_grid)-1)],mu[1,],type="l",col=color_palette[1],ylab="",xlab="")
for(i in seq(1,n,1)){
  lines(omega_grid[1:(length(omega_grid)-1)],mu[i,],type="l",col=color_palette[i])
}

# Write to CSV
write.csv(mu, "histograms.csv", row.names = FALSE)
write.csv(omega_grid[1:(length(omega_grid)-1)], "omega_grid.csv", row.names = FALSE)
