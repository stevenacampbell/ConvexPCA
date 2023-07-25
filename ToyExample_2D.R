#Reset Environment
rm(list=ls())

# Import source files
sourceCpp("HelperFunctions.cpp")
source("CPCAFunctions.R")

# Set seed
set.seed(123)

### Create data set ###

# Sample uniformly elements of the rectangle with sides [a,b] and [c,d]
sample_square <- function(a, b, c, d, num_samples) {
  x <- runif(num_samples, a, b)
  y <- runif(num_samples, c, d)
  data.frame(x = x, y = y)
}

# Sample data in a strip
num_samples <- 1000
w_1 <- 0.05; w_2 <- 0.25 # interval for x-axis
h_1 <- 0; h_2 <- 3.5 # interval for y-axis
dw<-w_2-w_1 # strip width

samples <- sample_square(w_1, w_2, h_1, h_2, num_samples)
X_init<-as.matrix(t(samples)) # initial data

# Create convex constraint set defined as the region to the right of the 
# vertical line x=0 and above the line y=m*x.
 
# Specify constraints
m_constr<-4
A<-rbind(c(1,0),c(-m_constr,1))
b_vec<-rep(0,nrow(A))

# Filter strip to only include data in the constraint set
id<-0
epsilon <- 0
for(i in seq(1,num_samples,1)){
  test<- A %*% X_init[,i]
  if(all(test>=epsilon)){
    if(id==0){
      id<-1
      X<-X_init[,i]
    }else{
      X<-cbind(X,X_init[,i])
    }
  }
}

# Add in adversarial points on the boundary and 
# near the lower left corner of the strip.

n_pts<-25 # number of adversarial points
delta<- 0.05*dw # perturbation parameter
perturbations <- delta * runif(n_pts) # sample random perturbations

# create extra adversarial points
extra_pts<-matrix(0,nrow=2,ncol=n_pts) 
extra_pts[1,]<-w_1+perturbations
extra_pts[2,]<-(w_1+perturbations)*m_constr

# add points to data set
X<-cbind(X,extra_pts)
x0<-as.vector(rowMeans(X)) # reference element

### Solve CPCA Problem ###

npcs<-1 # number of convex PCs
nthreads<-1 # number of threads to recruit
h<-0.000001 # central difference step size

CPCs<-CPCA(npcs,X,x0,A,b_vec,h,nthreads)
CPCs
PCs<-prcomp(t(X-x0),center=FALSE,scale=FALSE)
P<-PCs$rotation[,1] # allocate first PC
P

### Compute explained variation ###

TV<-sum((X-x0)^2) # Total variation of data

# Explained % variation of Euclidean PC
pts<-matrix(0,ncol=ncol(X),nrow=1)
for(i in seq(1,ncol(X),1)){
  pts[i]<-(X[,i]-x0) %*% P
}
per_var_pc1<-sum(pts^2)/TV

# Explained % variation of Convex PC
ts<-intersection_points(x0,CPCs,A,b_vec) # intersection points of PC with bdy
pts_cpc<-matrix(0,ncol=ncol(X),nrow=1)
for(i in seq(1,ncol(X),1)){
  pts_cpc[i]<-(X[,i]-x0) %*% CPCs
  test<- A %*% (x0+pts_cpc[i]*CPCs)
  if(!all(test>=0)){
    t_<-pts_cpc[i]
    if(abs(t_-ts[1])<=abs(t_-ts[2])){
      pts_cpc[i]<-ts[1]
    }else{
      pts_cpc[i]<-ts[2]
    }
  }
}
per_var_cpc1<-sum(pts_cpc^2)/TV

# Print results and observe CPC explained variation <= PC explained variation
print(per_var_pc1)
print(per_var_cpc1)

### Create Illustrative Figure ###

z<-seq(0,w_2+dw,0.001)
plot(seq(0,w_2+0.75*dw,(w_2+0.75*dw)/100), seq(0,h_2+0.05*h_2,(h_2+0.05*h_2)/100),
     xaxt = "n",yaxt = "n", xlab="", ylab="", type = "n")

# create a sequence of x values
x_seq = seq(0, 1000, length.out = 101)

# calculate corresponding y values
y_seq = m_constr * x_seq

# draw and fill the convex set
polygon(c(0, x_seq, max(x_seq)), c(max(y_seq), y_seq, max(y_seq)), col = "gray95")
box()

# draw lines
lines(x_seq, y_seq, lwd = 3, col="black")
lines(rep(0,length(y_seq)),seq(0,h_2+h_2*0.1,(h_2+h_2*0.1)/100), lwd = 3, col="black")

# Plot points and principal components
points(X[1,],X[2,],pch=0,cex=0.75,col="grey50")
m<-(x0[2]+CPCs[2]-x0[2])/(x0[1]+CPCs[1]-x0[1])
b<-x0[2]-m*x0[1]
lines(z[75:length(z)],m*z[75:length(z)]+b,col="blue",lwd=3)
m_<-(x0[2]+P[2]-x0[2])/(x0[1]+P[1]-x0[1])
b_<-x0[2]-m_*x0[1]
lines(z,m_*z+b_,col="purple",lwd=3,lty=2)
points(x0[1],x0[2],col="red",cex=1,pch=15)
lines(x_seq, y_seq, lwd = 4, col="black")
lines(rep(0,length(y_seq)),seq(0,h_2+0.1*h_2,(h_2+0.1*h_2)/100), lwd = 4, col="black")
legend("bottomright", 
       legend=c("Convex PC","Euclidean PC","Convex Set Bdy","Data Points", "Reference Point"), 
       col=c("blue","purple","black","grey50","red"), lty=c(1,2,1,NA,NA), 
       pch=c(NA,NA,NA,0,15), lwd=c(2,2,2,1,10), merge=FALSE, cex=0.75)


