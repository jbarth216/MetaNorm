#' Vector version of summing up every thing except for the ith element
#'
#' @param i The index to omit
#' @param vec The vector
#' @return The sum
#' @export
sum_except <- Vectorize(function(i,vec){sum(vec[-i])},vectorize.args="i")


#' Sum up every thing except for the ith element
#'
#' @param i The index to omit
#' @param vec The vector
#' @return The sum
#' @export
sum_except2 <- function(i,vec){sum(vec[-i])}


#' X
#'
#' @return A ratio
#' @export
X_thing <- function(i,j,X){
  X5 <- X[5]
  X6 <- X[6]
  ((X[i]-X6)*(X[j]-X6) + (X[i]-X5)*(X[j]-X5))/((X5-X6)^2)
}


#' Update a
#'
#' @return A posterior sample of a
#' @export
update_norm_a <- function(Y,X_pos,Y_neg,Y_hk,Y_reg,
                     cc,b,d_pos,d_neg,d_hk,d_reg,mu_a,phi,kappa_hk,kappa_reg,
                     sig2_a,sig2_e,sig2_n,P,N,H,R,n_samples){
  sigma2 <-  1/(1/sig2_a + (P+H+R)/sig2_e+N/sig2_n)
  A1 <- colSums(Y_neg - cc*matrix(b,N,n_samples,byrow=T) - matrix(d_neg,N,n_samples))
  A2 <- colSums(Y-matrix(X_pos[,2],P,1)%*%matrix(b,1,n_samples) -
                  matrix(d_pos,P,n_samples))
  A3 <- colSums(Y_hk-(matrix(phi,H,n_samples,byrow=T)+kappa_hk)*matrix(b,H,n_samples,byrow=T)
                -matrix(d_hk,H,n_samples))
  A4 <- colSums(Y_reg-(matrix(phi,R,n_samples,byrow=T)+kappa_reg)*matrix(b,R,n_samples,byrow=T)
                -matrix(d_reg,R,n_samples))
  mu <- sigma2*(mu_a/sig2_a + (A2+A3+A4)/sig2_e + A1/sig2_n)
  rnorm(n_samples,mu,sqrt(sigma2))
} ##done


#' Update b
#'
#' @return A posterior sample of b
#' @export
update_norm_b <- function(Y,X_pos,Y_neg,Y_hk,Y_reg,
                     cc,a,d_pos,d_neg,d_hk,d_reg,mu_b,phi,kappa_hk,kappa_reg,
                     sig2_b,sig2_e,sig2_n,P,N,H,R,n_samples){
  phi_hk <- matrix(phi,H,n_samples,byrow=T)+kappa_hk
  phi_reg <- matrix(phi,R,n_samples,byrow=T)+kappa_reg
  sigma2 <-1/(1/sig2_b +
                (sum(X_pos[,2]^2)+colSums(phi_hk^2)+colSums(phi_reg^2))/sig2_e
              +N*(cc^2)/sig2_n)
  B1 <- matrix(cc,1,N)%*%(Y_neg-matrix(a,N,n_samples,byrow=T) - matrix(d_neg,N,n_samples))
  B2 <- matrix(X_pos[,2],1,P)%*%(Y-matrix(a,P,n_samples,byrow=T) -
                                   matrix(d_pos,P,n_samples))
  B3 <- colSums(phi_hk*(Y_hk-matrix(a,H,n_samples,byrow=T)-matrix(d_hk,H,n_samples)))
  B4 <- colSums(phi_reg*(Y_reg-matrix(a,R,n_samples,byrow=T)-matrix(d_reg,R,n_samples)))
  mu <- sigma2*(mu_b/sig2_b + (B2+B3+B4)/sig2_e + B1/sig2_n)
  rnorm(n_samples,mu,sqrt(sigma2))
}


#' Update c
#'
#' @return A posterior sample of c
#' @export
update_norm_cc <- function(Y_neg,a,b,d_neg,sig2_n,L=-10,U=10,N,n_samples){
  sigma2 <- (N*sum(b^2)/sig2_n)^(-1)
  N1 <- sum(matrix(b,1,n_samples)%*%t(Y_neg - matrix(a,N,n_samples,byrow=T) - matrix(d_neg,N,n_samples)))
  mu <- sigma2*N1/sig2_n
  truncnorm::rtruncnorm(1,L,U,mu,sqrt(sigma2))
}


#' Update kappa
#'
#' @return A posterior sample of kappa
#' @export
update_norm_kappa <- function(Y_dat,a,b,phi,d,lambda,
                         sig2_e,sig2_kappa,n_genes,n_samples){
  sigma2 <- 1/(matrix(b^2,n_genes,n_samples,byrow=T)/sig2_e +
                 matrix(1/sig2_kappa,n_genes,n_samples))
  K1 <- ((Y_dat - matrix(a,n_genes,n_samples,byrow=T) - matrix(b*phi,n_genes,n_samples,byrow=T) -
            matrix(d,n_genes,n_samples))*matrix(b,n_genes,n_samples,byrow=T))/sig2_e
  K2 <- matrix(lambda,n_genes,n_samples)/sig2_kappa
  mu <- (K1 + K2)*sigma2
  matrix(rnorm(n_samples*n_genes,mu,sqrt(sigma2)),n_genes,n_samples)
}


#' Update d^-_n
#'
#' @return A posterior sample of d^-_n
#' @export
update_norm_dneg2 <- function(Y_neg,cc,a,b,d_neg,sig2_dn,sig2_n,n_samples,N,j){
  vj <- sum_except2(j,d_neg[1:(N-1)])
  sigma2 <- 1/(2/sig2_dn + 2*n_samples/sig2_n)
  D1<-sum(Y_neg[j,] - a - cc*b)
  D2 <- sum(Y_neg[N,] - a - cc*b + vj)
  mu <- sigma2*((D1-D2)/sig2_n - vj/sig2_dn)
  rnorm(1,mu,sqrt(sigma2))
}


#' Update d^+_n
#'
#' @return A posterior sample of d^+_n
#' @export
update_norm_dpos3 <- function(Y,X_pos,a,b,sig2_d,sig2_e,P,n_samples,Sig_uns,D5,D6){
  Ytil <- Y-matrix(a,P,n_samples,byrow=T)- outer(X_pos[,2],b)
  scaler <- (1/sig2_d + n_samples/sig2_e)^(-1)
  Sig <- scaler*Sig_uns
  Phi <- apply((Ytil[1:4,] - D6%*%Ytil[5,] + D5%*%Ytil[6,]),1,sum)/sig2_e
  mu <- Sig%*%Phi
  d<-mvtnorm::rmvnorm(1,mu,Sig)
  d[5]<- -t(D6)%*%t(d)
  d[6]<- -sum(d)
  d
}


#' Update d^*_n
#'
#' @return A posterior sample of d^*_n
#' @export
update_norm_dhkreg <- function(Y,a,b,phi,kappa,sig2_e,sig2_d,n_genes,n_samples){
  sigma2 <- 1/(n_samples/sig2_e + 1/sig2_d)
  phi_kappa <- matrix(phi,n_genes,n_samples,byrow=T)+kappa
  mu <- sigma2*rowSums(Y-matrix(a,n_genes,n_samples,byrow=T)
                       -matrix(b,n_genes,n_samples,byrow=T)*phi_kappa)/sig2_e
  rnorm(n_genes,mu,sqrt(sigma2))
}


#' Update mu
#'
#' @return A posterior sample of mu
#' @export
update_norm_mu <- function(z,sig2_z,prior_mean,prior_var,L=-10,U=10,n_samples){
  sigma2 <- 1/(n_samples/sig2_z + 1/prior_var)
  mu <- sigma2*(sum(z)/sig2_z + prior_mean/prior_var)
  truncnorm::rtruncnorm(1,L,U,mu,sqrt(sigma2))
}


#' Update sigma^2_a, sigma^2_b, and sigma^2_d
#'
#' @return A posterior sample of sigma^2_a, sigma^2_b, and sigma^2_d
#' @export
update_norm_sig2_abd <- function(z,mu_z,n,eps=.01){
  alpha <- eps+n/2
  beta <- eps + sum((z-mu_z)^2)/2
  1/rgamma(1,alpha,beta)
}


#' Update sigma^2_e-
#'
#' @return A posterior sample of sigma^2_e-
#' @export
update_norm_sig2_n <- function(Y_neg,a,b,cc,d_neg,n_genes,n_samples,eps=.01){
  alpha <- eps + (n_genes*n_samples)/2
  beta <- eps + sum((Y_neg-matrix(a,n_genes,n_samples,byrow=T)-
                       matrix(b,n_genes,n_samples,byrow=T)*cc-
                       matrix(d_neg,n_genes,n_samples))^2)/2
  1/rgamma(1,alpha,beta)
}


#' Update sigma^2_e
#'
#' @return A posterior sample of sigma^2_e
#' @export
update_norm_sig2_e <- function(Y_pos,X_pos,Y_hk,Y_reg,a,b,phi,kappa_hk,kappa_reg,
                          d_pos,d_hk,d_reg,P,H,R,n_samples,eps=.01){
  phi_hk <- matrix(phi,H,n_samples,byrow=T)+kappa_hk
  phi_reg <- matrix(phi,R,n_samples,byrow=T)+kappa_reg
  alpha <- eps + ((P+H+R)*n_samples)/2
  beta <- eps + (sum((Y_pos-matrix(a,P,n_samples,byrow=T)-
                        matrix(X_pos[,2],P,1)%*%matrix(b,1,n_samples) -
                        matrix(d_pos,P,n_samples))^2)
                 + sum((Y_hk-matrix(a,H,n_samples,byrow=T) -
                          matrix(b,H,n_samples,byrow=T)*phi_hk -
                          matrix(d_hk,H,n_samples))^2)
                 + sum((Y_reg-matrix(a,R,n_samples,byrow=T) -
                          matrix(b,R,n_samples,byrow=T)*phi_reg -
                          matrix(d_reg,R,n_samples))^2) )/2
  1/rgamma(1,alpha,beta)
}


#' Update sigma^2_kappa
#'
#' @return A posterior sample of sigma^2_kappa
#' @export
update_norm_sig2_kappa <- function(kappa,lambda,n_genes,n_samples,eps=.01){
  alpha <- eps + (n_genes*n_samples)/2
  beta <- eps + sum((kappa - matrix(lambda,n_genes,n_samples))^2)/2
  1/rgamma(1,alpha,beta)
}

