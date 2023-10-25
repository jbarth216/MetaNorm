library(MASS)
library(truncnorm)


#' Update parameter a_{ik} and b_{ik}
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter \theta_{ik} = (a_{ik}, b_{ik})
#' @param i The patient index
#' @param k The study index
#' @param indx The prob index
#' @param m The current posterior sample of the mean
#' @param Sig The current posterior sample of the covariance
#' @param sigjk The current posterior sample of the variances of the residuals
#' @param sjk The current posterior sample of s_{jk}
#' @return A sample of the \theta parameter
#' @export
update_ab<-Vectorize(function(i,k,indx,m,Sig,sigjk,sjk){
  m <- m[,k]
  Sig <- Sig[[k]]
  sigjk <- 1/sigjk[,k]
  Y <- Yik_list[[indx]]$Y
  s <- sjk[,k]
  # We first compute the updated covariance
  Sig1 <- Sig1_cpp(Sig,X,diag(sigjk))
  # Using the updated covariance, we update the mean
  mu <- mu_cpp(Sig1,Sig,X,diag(sigjk),m,Y,s)
  rMvNorm(1,mu,Sig1)
},vectorize.args=c("i","k","indx"))


#' Update parameter Sigma_k
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter Sigma_k
#' @param k The study index
#' @param ab1 The current posterior samples of the intercept and the slope
#' @param m The current posterior sample of the mean vector
#' @return A sample of the Sigma parameter
#' @export
update_Sig <- function(k,ab1,m){
  mu_rep<-matrix(rep(m[,k],nk[k]),nk[k],2,byrow=T)
  Phi<-as.matrix(t(ab1[ab1$k==k,3:4]-mu_rep))%*%as.matrix((ab1[ab1$k==k,3:4]-mu_rep))
  RateMat<-inverse(matrix(I2+Phi,2,2))
  inverse(matrix(rWishart(1,3+nk[k],RateMat),2,2))
}


#' Update parameter m
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter m
#' @param k The study index
#' @param ab The current posterior sample of the intercept and the slope
#' @param sig_alpha The current posterior sample of sigma^2_alpha
#' @param sig_beta The current posterior sample of sigma^2_beta
#' @param Sig The current posterior sample of the covarance
#' @param mu The current posterior sample of the mu parameter
#' @return A sample of the m parameter
#' @export
update_m <- Vectorize(function(k,ab,sig_alpha,sig_beta,Sig,mu){
  theta <- t(as.matrix(ab[ab$k==k,3:4]))
  Sig_m <- matrix(c(sig_alpha,0,0,sig_beta),2,2)
  Sig_m_inv <- solve(Sig_m)
  Sig_inv <- solve(Sig[[k]])
  Var <- solve(Sig_m_inv+nk[k]*Sig_inv)
  Mean <- Var%*%(Sig_m_inv%*%mu+apply(Sig_inv%*%theta,1,sum))
  mvrnorm(1,Mean,Var)
},vectorize.args="k")


#' Update parameter mu
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter mu. As the form of the updating rules
#' are identical for mu_alpha and mu_beta, we use the same
#' function for both of them
#' @param L The lower bound
#' @param U The upper bound
#' @param alpha The current posterior sample of alpha/beta
#' @param sig_alpha The current posterior sample of sigma^2_alpha/sigma^2_beta
#' @param K The number of studies
#' @return A sample of the mu_alpha/mu_beta parameter
#' @export
update_mu <- function(L,U,alpha,sig_alpha,K=13){
  if (L >= U){stop("Lower bound cannot be larger than Upper bound")}
  Mean <- (1/K)*sum(alpha)
  Sd <- ((1/K)*sig_alpha)^.5
  rtruncnorm(1,L,U,Mean,Sd)
}


#' Update parameter sigma_alpha/sigma_beta
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter sigma_alpha/sigma_beta As the form of the updating rules
#' are identical for sigma_alpha and sigma_beta, we use the same
#' function for both of them
#' @param eps The IG parameter
#' @param alpha The current posterior sample of alpha/beta
#' @param mu The current posterior sample of mu_alpha/mu_beta
#' @param K The number of studies
#' @return A sample of the sigma_alpha/sigma_beta parameter
#' @export
update_sig <- function(eps,alpha,mu,K=13){
  a <- eps + (K/2)
  b <- eps + .5*sum((alpha-mu)^2)
  1/rgamma(1,a,b)
}


#' Update parameter s_jk
#'
#' This function performs one step of the Gibbs sampling
#' of the first 4 components of the parameter s_jk
#' @param j The prob index
#' @param k The study index
#' @param tj The current posterior sample of t_j
#' @param tau2j The current posterior sample of tau^2_j
#' @param sigjk The current posterior sample of sigma^2_jk
#' @param ab The current posterior sample of theta
#' @param sjk The current posterior sample of sjk
#' @return A sample of the sjk parameter
#' @export
update_sjk <- Vectorize(function(j,k,tj,tau2j,sigjk,ab,sjk){
  Sj <- sum(sjk[c(-j,-5,-6),k])
  XSj <- sum(t(X[c(-j,-5,-6),2])%*%sjk[c(-j,-5,-6),k])
  uj <- (X[6,2]-X[j,2])/(X[5,2]-X[6,2])
  mm <- (X[6,2]*Sj - XSj)/(X[5,2]-X[6,2])
  Var <- ((1/tau2j[j])+(uj^2)/tau2j[5]+((1+uj)^2)/tau2j[6]+
            (nk[k]/sigjk[j,k])+(nk[k]*uj^2)/sigjk[5,k]+
            (nk[k]*(1+uj)^2)/sigjk[6,k])^(-1)
  Yj<-Yjk_list[[k]][[j]]$Y
  Y5<-Yjk_list[[k]][[5]]$Y
  Y6<-Yjk_list[[k]][[6]]$Y
  abk <- ab[ab$k==k,3:4]
  Mean <- Var*((tj[j]/tau2j[j])+((tj[5]-mm)*uj/tau2j[5])-
                 (tj[6]+Sj+mm)*(1+uj)/tau2j[6] +
                 (1/sigjk[j,k])*sum(Yj-X[j,]%*%t(abk))+
                 (1/sigjk[5,k])*sum(Y5-X[5,]%*%t(abk)-mm)*uj -
                 (1/sigjk[6,k])*sum(Y6-X[6,]%*%t(abk)+Sj+mm)*(1+uj))
  rnorm(1,Mean,sqrt(Var))
},vectorize.args=c("j","k"))


#' Update parameter s_5k
#'
#' This function performs one step of the Gibbs sampling
#' of the 5th components of the parameter s_jk
#' @param k The study index
#' @param sjk The current posterior sample of sjk
#' @return A sample of the s5k parameter
#' @export
update_s5k <- Vectorize(function(k,sjk){
  s <- sum(sjk[1:4,k])
  Xs <- sum(t(X[1:4,2])%*%sjk[1:4,k])
  (X[6,2]*s - Xs)/(X[5,2]-X[6,2])
},vectorize.args="k")


#' Update parameter tj
#'
#' This function performs one step of the Gibbs sampling
#' of the first 4 components of the parameter tj
#' @param j The prob index
#' @param L The lower bound
#' @param U The upper bound
#' @param sjk The current posterior sample of sjk
#' @param tau2j The current posterior sample of tau^2_j
#' @param tj The current posterior sample of t_j
#' @param K The number of studies
#' @return A sample of the tj parameter
#' @export
update_tj <- Vectorize(function(j,L,U,sjk,tau2j,tj,K=13){
  if (L >= U){stop("Lower bound cannot be larger than Upper bound")}
  Tj <- sum(tj[c(-j,-5,-6)])
  XTj <- sum(t(X[c(-j,-5,-6),2])%*%tj[c(-j,-5,-6)])
  uj <- (X[6,2]-X[j,2])/(X[5,2]-X[6,2])
  mm <- (X[6,2]*Tj - XTj)/(X[5,2]-X[6,2])
  Var <- (1/K)*((1/tau2j[j])+(uj^2)/tau2j[5] + ((1+uj)^2)/tau2j[6])^(-1)
  Mean <-  Var*sum(sjk[j,]/tau2j[j]+(uj/tau2j[5])*(sjk[5,]-mm) -
                     ((1+uj)/tau2j[6])*(sjk[6,]+Tj+mm))
  rtruncnorm(1,L,U,Mean,sqrt(Var))
},vectorize.args="j")


#' Update parameter tau2j
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter tau2j
#' @param j The prob index
#' @param eps The IG paramter
#' @param sjk The current posterior sample of sjk
#' @param tj The current posterior sample of t_j
#' @param K The number of studies
#' @return A sample of the tau2j parameter
#' @export
update_tau2j <- Vectorize(function(j,eps,sjk,tj,K=13){
  a <- K/2 + eps
  b <- eps + .5*sum((sjk[j,]-tj[j])^2)
  1/rgamma(1,a,b)
},vectorize.args="j")


#' Update parameter sigma^2_jk
#'
#' This function performs one step of the Gibbs sampling
#' of the parameter sigma^2_jk
#' @param j The prob index
#' @param k The study index
#' @param eps The IG paramter
#' @param ab The current posterior sample of theta
#' @param sjk The current posterior sample of sjk
#' @return A sample of the sigma^2_jk parameter
#' @export
update_sigjk <- Vectorize(function(j,k,eps,ab,sjk){
  Y <- Yjk_list[[k]][[j]]$Y
  a <- eps + nk[k]/2
  b <- eps + .5*sum((Y-X[j,]%*%t(as.matrix(ab[ab$k==k,3:4]))-sjk[j,k])^2)
  1/rgamma(1,a,b)
},vectorize.args=c("j","k"))


