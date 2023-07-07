#include <RcppArmadillo.h>
#include <omp.h>  // Include OpenMP library
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec undo_spherical_coords(arma::vec theta){
  int len = theta.n_elem;
  arma::vec P(len + 1, arma::fill::zeros);  // create a vector of size len+1, and fill it with zeros
  if (len == 1){
    P[0] = cos(theta[0]);
    P[1] = sin(theta[0]);
  } else {
    arma::vec stheta = sin(theta);
    arma::vec sthetacp = cumprod(stheta);
    arma::vec ctheta = cos(theta);
    arma::vec fsin(len + 1, arma::fill::zeros);
    arma::vec fcos(len + 1, arma::fill::zeros);
    
    fsin[0] = 1;
    fcos[len] = 1;
    
    for (int i = 0; i < len + 1; ++i){
      if (i < len) {
        fcos[i] = ctheta[i];
        fsin[i + 1] = sthetacp[i];
      }
      P[i] = fcos[i] * fsin[i];
    }
  }
  return P;
}


// [[Rcpp::export]]
arma::vec spherical_coords(arma::vec P) {
  int d = P.n_elem;
  arma::vec theta(d - 1);
  
  if (d == 2) {
    theta[0] = acos(P[0]);
  } else {
    arma::vec P2 = square(P);
    arma::vec P2r = flipud(P2);
    arma::vec CS = cumsum(P2r);
    arma::vec CSr = flipud(CS);
    arma::vec CSr_ = CSr.subvec(0, d - 3);
    arma::vec P_ = P.subvec(0, d - 3);
    
    for (int i = 0; i < d - 2; i++) {
      theta[i] = acos(P_[i] / sqrt(CSr_[i]));
    }
    
    if (P[d - 1] < 0) {
      theta[d - 2] = 2 * arma::datum::pi - acos(P[d - 2] / sqrt(CSr[d - 2]));
    } else {
      theta[d - 2] = acos(P[d - 2] / sqrt(CSr[d - 2]));
    }
  }
  
  return theta;
}


// [[Rcpp::export]]
arma::mat project_orthogonal_complement(arma::mat X, arma::mat B) {
  // Calculate the projection of X onto B
  arma::mat projection_B = B * B.t() * X;
  
  // Subtract the projection from X to get the projection onto the orthogonal complement
  arma::mat projection_B_orthogonal = X - projection_B;
  
  return projection_B_orthogonal;
}


// [[Rcpp::export]]
arma::vec ls_fit(arma::mat B, arma::vec p_) {
  arma::vec p = solve(B, p_);
  return p;
}

// [[Rcpp::export]]
arma::vec intersection_points(arma::vec x0, arma::vec p, arma::mat A, arma::vec b) {
  int dim = b.n_elem;
  arma::vec tvals = arma::vec(2).fill(arma::datum::inf);
  arma::vec num = b - A * x0;
  arma::vec den = A * p;
  
  for(int i = 0; i < dim; ++i){
    if(std::abs(den[i]) > std::numeric_limits<double>::epsilon()){
      double a = num[i] / den[i];
      arma::vec x = x0 + a * p;
      arma::vec constraint_val = A * x - b;
      
      // To implement "round", we're going to do:
      constraint_val = arma::round(constraint_val * 1e10) / 1e10;
      if(arma::all(constraint_val >= 0)){
        if(tvals[0] == arma::datum::inf){
          tvals[0] = a;
        } else {
          if(std::abs(a - tvals[0]) > std::numeric_limits<double>::epsilon()){
            tvals[1] = a;
            break;
          }
        }
      }
    }
  }
  return tvals;
}

// [[Rcpp::export]]
double restricted_projection_coord(arma::vec x, arma::vec x0, arma::vec p, arma::mat A, arma::vec b, arma::vec tvals) {
  double a = dot(x - x0, p);
  arma::vec x_ = x0 + a * p;
  arma::vec constraint_val = A * x_ - b;
  
  // To implement "round", we're going to do:
  constraint_val = arma::round(constraint_val * 1e10) / 1e10;
  double t_sol;
  
  if(arma::all(constraint_val>=0)){
    t_sol = a;
  } else {
    double t1 = tvals[0];
    double t2 = tvals[1];
    if(std::abs(a-t1) < std::abs(a-t2)){
      t_sol = t1;
    } else {
      t_sol = t2;
    }
  }
  return t_sol;
}

// [[Rcpp::export]]
double eval_CPCA_obj(arma::vec theta, arma::mat B, arma::vec x0, arma::mat A, arma::vec b, arma::mat X) {
  arma::vec omega = undo_spherical_coords(theta);
  arma::vec p = B * omega;
  arma::vec tvals = intersection_points(x0, p, A, b);
  int nobs = X.n_cols;
  double val = 0;
  
  for(int i = 0; i < nobs; ++i){
    arma::vec x = X.col(i);
    double t = restricted_projection_coord(x, x0, p, A, b, tvals);
    arma::vec sol = x0 + t * p;
    val += arma::sum(arma::pow(x - sol, 2));
  }
  
  val /= nobs;
  
  return(val);
}

// [[Rcpp::export]]
arma::vec eval_grad_CPCA_obj(arma::vec theta, double h, arma::mat B, arma::vec x0, arma::mat A, arma::vec b, arma::mat X) {
  int d = theta.n_elem;
  arma::vec grad(d, arma::fill::zeros);
  
  for(int i = 0; i < d; ++i){
    arma::vec theta_p = theta;
    arma::vec theta_m = theta;
    theta_p[i] += h;
    theta_m[i] -= h;
    grad[i] = (eval_CPCA_obj(theta_p, B, x0, A, b, X) - eval_CPCA_obj(theta_m, B, x0, A, b, X)) / (2.0 * h);
  }
  
  return grad;
}


// [[Rcpp::export]]
arma::vec eval_grad_CPCA_obj_parallel(arma::vec theta, double h, arma::mat B, arma::vec x0, arma::mat A, arma::vec b, arma::mat X, int num_threads) {
  
  // Set the number of threads
  omp_set_num_threads(num_threads);
  
  int d = theta.n_elem;
  arma::vec grad(d);
  
  #pragma omp parallel for
  for(int i = 0; i < d; ++i){
    arma::vec theta_p = theta;
    arma::vec theta_m = theta;
    theta_p[i] += h;
    theta_m[i] -= h;
    grad[i] = (eval_CPCA_obj(theta_p, B, x0, A, b, X) - eval_CPCA_obj(theta_m, B, x0, A, b, X)) / (2.0 * h);
  }
  
  return grad;
}