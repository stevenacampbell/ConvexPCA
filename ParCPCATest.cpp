#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

std::vector<double> MatrixMult_v(RcppParallel::RVector<double> Ad, RcppParallel::RVector<double> Ad_, std::vector<double> x, int dim){
  std::vector<double> v(dim);
  for(int i=0; i<dim; i++){
    if(i==0){
      v[i]=Ad[i]*x[i];
    }else{
      v[i]=Ad_[i-1]*x[i-1]+Ad[i]*x[i];
    }
  }
  return v;
}


double gett0_v2(RcppParallel::RVector<double> X0, std::vector<double> Xd, std::vector<double> P, int dim){
  double t0=0;
  for(int i = 0; i<dim; i++){
    t0 = t0+(Xd[i]-X0[i])*P[i];
  }
  return t0;
}

std::vector<double> GetCol(RcppParallel::RMatrix<double> X ,int i, int d){
  RMatrix<double>::Column Xcol = X.column(i);
  std::vector<double> Xd(d);
  for(int j=0; j<d; ++j){
    Xd[j]=Xcol[j];
  }
  return Xd;
}

std::vector<double> GetRow(RcppParallel::RMatrix<double> X ,int i, int d){
  RMatrix<double>::Row Xrow = X.row(i);
  std::vector<double> Xd(d);
  for(int j=0; j<d; ++j){
    Xd[j]=Xrow[j];
  }
  return Xd;
}

std::vector<double> MatrixMult_B(RcppParallel::RMatrix<double> B, std::vector<double> x, int dim){
  std::vector<double> v(dim);
  std::vector<double> z(dim);
  std::vector<double> b(dim);
  for(int i=0;i<dim;i++){
    b = GetRow(B,i,dim);
    z[i]=0;
    for(int j=0;j<dim;j++){
      z[i]=z[i]+b[j]*x[j];
    }
  }
  return z;
}

std::vector<double> MatrixMult_B_v1(RcppParallel::RMatrix<double> B, std::vector<double> x, int dim, int ncol){
  std::vector<double> z(dim);
  std::vector<double> b(ncol);
  for(int i=0;i<dim;i++){
    b = GetRow(B,i,ncol);
    z[i]=0;
    for(int j=0;j<ncol;j++){
      z[i]=z[i]+b[j]*x[j];
    }
  }
  return z;
}

double Projection_v2(RcppParallel::RVector<double> Ad, RcppParallel::RVector<double> Ad_, RcppParallel::RVector<double> X0, std::vector<double> Xd, std::vector<double> P, RcppParallel::RVector<double>  NUM, int dim, std::vector<double> tvals, std::vector<int> id){
  double t;
  double tsol = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  std::vector<double> z(dim);
  int test=0;
  double eps = (-1)*pow(10,-10);
  std::vector<double> x(dim);
  
  t = gett0_v2(X0,Xd,P,dim);
  for(int i=0;i<dim;i++){
    z[i]=X0[i]+t*P[i];
  }
  x = MatrixMult_v(Ad,Ad_,z,dim);
  //Check if all positive up to numerical error
  for(int i=0;i<dim;i++){
    if(x[i]<eps){
      test=1;
      break;
    }
  }
  if(test==0){
    tsol = t;
  }else{
    if(id[1]==-1){
      tsol=tvals[id[0]-1];
    }else{
      t1 = tvals[id[0]-1];
      t2 = tvals[id[1]-1];
      if(abs(t-t1)<abs(t-t2)){
        tsol = t1;
      }else{
        tsol = t2;
      }
    }
  }
  return tsol;
}


double ObjFuncEval_v2(RcppParallel::RVector<double> Ad, RcppParallel::RVector<double> Ad_, RcppParallel::RVector<double> X0, std::vector<double> P, RcppParallel::RVector<double>  NUM, RcppParallel::RMatrix<double> X, int nobs, double inobs, int d){
  double obj = 0.0;
  double t = 0.0;
  double sum = 0.0;
  double sol = 0.0;
  std::vector<double> Xd(d);
  std::vector<double> tvals(d);
  std::vector<int> id(2);
  id[0]=-1;
  id[1]=-1;
  int test1;
  int ind=0;
  double eps = (-1)*pow(10,-10);
  double epsp = pow(10,-10);
  std::vector<double> DEN(d);
  std::vector<double> x1(d);
  std::vector<double> v(d);
  
  
  for(int i=0;i<nobs;++i){
    Xd = GetCol(X,i,d);
    DEN = MatrixMult_v(Ad,Ad_,P,d);
    for (int i=0; i < d; ++i) {
      if(DEN[i]!=0){
        tvals[i]= -(NUM[i]/DEN[i]);
        for(int j=0; j<d; ++j){
          v[j] = X0[j]-(NUM[i]/DEN[i])*P[j];
        }
        x1=MatrixMult_v(Ad,Ad_,v,d);
        test1=0;
        for(int j=0; j<d; j++){
          if(x1[j]<eps){
            test1=1;
            break;
          }
        }
        if(test1==0){
          if(ind==0){
            id[ind]=i+1;
            ind=ind+1;
          }else{
            if(abs(tvals[i]-tvals[id[0]-1])>epsp){
              id[ind]=i+1;
              ind=ind+1;
            }
          }
        }
        if(ind>1){
          break;
        }
      }
    }
    t = Projection_v2(Ad, Ad_, X0, Xd, P, NUM, d, tvals, id);
    sum = 0.0;
    for(int j=0;j<d;++j){
      sol = X0[j]+t*P[j];
      sum = sum + pow(sol-Xd[j],2);
    }
    obj = obj+inobs*sum;
  }
  return obj;
}


std::vector<double> CP(std::vector<double> v, int len){
  std::vector<double> v_(len);
  for(int i=0; i<len; i++){
    if(i==0){
      v_[i]=v[i];
    }else{
      v_[i]=v[i]*v_[i-1];
    }
  }
  return v_;
}

std::vector<double> GetP(RcppParallel::RMatrix<double> B, std::vector<double> theta, int len){
  int l1 = len+1;
  std::vector<double> P(l1);
  std::vector<double> stheta(len);
  std::vector<double> sthetacp(len);
  std::vector<double> ctheta(len);
  std::vector<double> fsin(l1);
  std::vector<double> fcos(l1);
  if(len == 1){
    P[0]=cos(theta[0]);
    P[1]=sin(theta[0]);
  }else{
    for(int j=0;j<len;j++){
      stheta[j] = sin(theta[j]);
      ctheta[j] = cos(theta[j]);
    }
    sthetacp = CP(stheta, len);
    fsin[0]=1;
    fcos[len]=1;
    for(int i=0; i<len+1;++i){
      if(i<len){
        fcos[i]=ctheta[i];
        fsin[i+1]=sthetacp[i];
      }
      P[i]=fcos[i]*fsin[i];
    }
  }
  P=MatrixMult_B(B, P, l1);
  return P;
}

std::vector<double> GetPNew(RcppParallel::RMatrix<double> B, std::vector<double> theta, int dim, int len){
  int l1 = len+1;
  std::vector<double> P(l1);
  std::vector<double> stheta(len);
  std::vector<double> sthetacp(len);
  std::vector<double> ctheta(len);
  std::vector<double> fsin(l1);
  std::vector<double> fcos(l1);
  std::vector<double> P0(dim);
  if(len == 1){
    P[0]=cos(theta[0]);
    P[1]=sin(theta[0]);
  }else{
    for(int j=0;j<len;j++){
      stheta[j] = sin(theta[j]);
      ctheta[j] = cos(theta[j]);
    }
    sthetacp = CP(stheta, len);
    fsin[0]=1;
    fcos[len]=1;
    for(int i=0; i<len+1;++i){
      if(i<len){
        fcos[i]=ctheta[i];
        fsin[i+1]=sthetacp[i];
      }
      P[i]=fcos[i]*fsin[i];
    }
  }
  P0=MatrixMult_B_v1(B, P, dim, l1);
  return P0;
}

//Define Derivative Approximation
double DerApprox(RcppParallel::RMatrix<double> B, RcppParallel::RVector<double> Ad, RcppParallel::RVector<double> Ad_, RcppParallel::RVector<double> X0, RcppParallel::RVector<double> theta, RcppParallel::RVector<double> NUM, RcppParallel::RMatrix<double> X, int nobs, double inobs, int d, int len, double h, int i){
  std::vector<double> thetap(len);
  std::vector<double> thetam(len);
  std::vector<double> Pp(len+1);
  std::vector<double> Pm(len+1);
  double valp;
  double valm;
  
  for (int j = 0; j < len; j++){
    thetap[j]=theta[j];
    thetam[j]=theta[j];
    if(j==i){
      thetap[j]=thetap[j]+h;
      thetam[j]=thetam[j]-h;
    }
  }
  
  Pp = GetPNew(B, thetap, d, len);
  Pm = GetPNew(B, thetam, d, len);
  valp = ObjFuncEval_v2(Ad, Ad_, X0, Pp, NUM, X, nobs, inobs, d);
  valm = ObjFuncEval_v2(Ad, Ad_, X0, Pm, NUM, X, nobs, inobs, d);
  
  return (valp-valm)/(2*h);
}

//Compute Gradient in Parallel
struct ParGrad : public Worker {
  
  // input matrix to read from
  const RMatrix<double> B;
  const RVector<double> Ad; 
  const RVector<double> Ad_;
  const RVector<double> X0;
  const RVector<double> theta; 
  const RVector<double> NUM;
  const RMatrix<double> X;
  int nobs;
  double inobs;
  int d;
  int len;
  double h;
  
  // output matrix to write to
  RMatrix<double> outvec;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  ParGrad(const NumericMatrix B, const NumericVector Ad, const NumericVector Ad_, const NumericVector X0, const NumericVector theta, const NumericVector NUM, const NumericMatrix X, int nobs, double inobs, int d, int len, double h, NumericMatrix outvec)
    : B(B), Ad(Ad), Ad_(Ad_), X0(X0), theta(theta), NUM(NUM), X(X), nobs(nobs), inobs(inobs), d(d), len(len), h(h), outvec(outvec) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      
      // write to output matrix
      outvec(i,0) = DerApprox(B, Ad, Ad_, X0, theta, NUM, X, nobs, inobs, d, len, h, i);
    }
  }
};

// [[Rcpp::export]]
NumericMatrix GradientEval(NumericMatrix B, NumericMatrix Ad, NumericVector Ad_, NumericVector X0, NumericVector theta, NumericVector NUM, NumericMatrix X, int nobs, double inobs, int d, int len, double h) {
  
  // allocate the matrix we will return
  NumericMatrix outvec(len, 1);
  
  // create the worker
  ParGrad parGrad(B, Ad, Ad_, X0, theta, NUM, X, nobs, inobs, d, len, h, outvec);
  
  // call it with parallelFor
  parallelFor(0, len, parGrad);
  
  return outvec;
}

// [[Rcpp::export]]
double FunctionEval(NumericMatrix B, NumericVector Ad, NumericVector Ad_, NumericVector X0, NumericVector theta, NumericVector NUM, NumericMatrix X, int nobs, double inobs, int d, int len) {
  
  // allocate the matrix we will return
  double outval;
  std::vector<double> P(d);
  std::vector<double> theta0(len);
  
  for (int j = 0; j < len; j++){
    theta0[j]=theta[j];
  }
  
  P = GetPNew((RcppParallel::RMatrix<double>) B, theta0, d, len);
  outval = ObjFuncEval_v2((RcppParallel::RVector<double>) Ad, (RcppParallel::RVector<double>) Ad_, (RcppParallel::RVector<double>) X0, P, (RcppParallel::RVector<double>) NUM, (RcppParallel::RMatrix<double>) X, nobs, inobs, d);
  return outval;
}

NumericVector CP_v(NumericVector v, int len){
  NumericVector v_(len);
  for(int i=0; i<len; i++){
    if(i==0){
      v_[i]=v[i];
    }else{
      v_[i]=v[i]*v_[i-1];
    }
  }
  return v_;
}

// [[Rcpp::export]]
NumericVector GetP_v(NumericVector theta, int len){
  int l1 = len+1;
  NumericVector P(l1);
  NumericVector stheta(len);
  NumericVector sthetacp(len);
  NumericVector ctheta(len);
  NumericVector fsin(l1);
  NumericVector fcos(l1);
  if(len == 1){
    P[0]=cos(theta[0]);
    P[1]=sin(theta[0]);
  }else{
    stheta = sin(theta);
    sthetacp = CP_v(stheta, len);
    ctheta = cos(theta);
    fsin[0]=1;
    fcos[len]=1;
    for(int i=0; i<len+1;++i){
      if(i<len){
        fcos[i]=ctheta[i];
        fsin[i+1]=sthetacp[i];
      }
      P[i]=fcos[i]*fsin[i];
    }
  }
  return P;
}


