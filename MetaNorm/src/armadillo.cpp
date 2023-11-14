// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


//' Generate random variables from a multivariate Gaussian distribution
//'
//' This function uses the Cholesky decomposition to transform
//' isotropic Gaussian samples to samples from the desired Gaussian distribution
//' @param n The number of samples
//' @param mu The mean vector
//' @param sigma The variance-covariance matrix
//' @return A matrix of random draws from the specified Gaussian distribution
//' @export
// [[Rcpp::export]]
arma::mat rMvNorm(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


//' Compute matrix inverse
//'
//' This function computes the inverse of the supplied matrix
//' @param sig A matrix to be inversed
//' @return The inverse of the supplied matrix
//' @export
// [[Rcpp::export]]
arma::mat inverse(arma::mat sig){
  return inv(sig);
}


//' A helper function of the Gibbs sampler to the Sigma parameters of a_{ik} and b_{ik}
//'
//' This function computes the posterior covariance matrix
//' of the parameter theta_{ik} = (a_{ik} and b_{ik})
//' @param Sig The current posterior sample of the covariance matrix
//' @param X A matrix of X_{ik}
//' @param sigjk A diagonal matrix of the posterior samples of sigma^2_{jk}
//' @return An updated Sigma matrix
//' @export
// [[Rcpp::export]]
arma::mat Sig1_cpp(arma::mat Sig, arma::mat X, arma::mat sigjk){
  return inv(inv(Sig)+X.t()*sigjk*X);
}


//' A helper function of the Gibbs sampler to mu the parameters a_{ik} and b_{ik}
//'
//' This function computes the posterior mean vector
//' of the parameter theta_{ik} = (a_{ik} and b_{ik})
//' @param Sig1 The updated posterior covariance matrix
//' @param Sig The current posterior sample of the covariance matrix
//' @param X X A matrix of X_{ik}
//' @param sigjk A diagonal matrix of the posterior samples of sigma^2_{jk}
//' @param m The current posterior sample of the mean vector
//' @param Y A vector of Y_{ik}
//' @param s A vector of s_k
//' @return An updated mu vector
//' @export
// [[Rcpp::export]]
arma::mat mu_cpp(arma::mat Sig1, arma::mat Sig,
                 arma::mat X, arma::mat sigjk, arma::vec m,
                 arma::vec Y, arma::vec s){
  return Sig1*(inv(Sig)*m+X.t()*sigjk*(Y-s));
}
