// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat rMvNorm(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat inverse(arma::mat sig){
  return inv(sig);
}

// [[Rcpp::export]]
arma::mat Sig1_cpp(arma::mat Sig, arma::mat X, arma::mat sigjk){
  return inv(inv(Sig)+X.t()*sigjk*X);
}

// [[Rcpp::export]]
arma::mat mu_cpp(arma::mat Sig1, arma::mat Sig,
                 arma::mat X, arma::mat sigjk, arma::vec m,
                 arma::vec Y, arma::vec s){
  return Sig1*(inv(Sig)*m+X.t()*sigjk*(Y-s));
}
