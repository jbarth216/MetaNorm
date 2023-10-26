#' Perform meta analysis
#'
#' This function performs meta analysis via Gibbs sampling
#' @param ds A curated data.frame
#' @param coeffs2 A data.frame of estimated coefficients
#' @param M Number of draws
#' @param n_keep Number of samples to keep. If <= 0, all samples will be kept.
#' @return A data.frame with posterior draws
#' @export
meta_analysis <- function(ds, coeffs2, M=5000, n_keep=-1)
{
  if((n_keep <= 0) | (n_keep > M))
  {
    n_keep = M
  }
  ##Some items that won't change##
  #####################################
  # AS THE GIBBS NEED THOSE VARIABLE
  # I AM USING <<-
  # THIS IS SUPER MESSY, SHOULD BE CHANGED
  #####################################
  I2 <<- matrix(c(1,0,0,1),2,2)
  X <<- matrix(c(rep(1,6), log(c(128,32,8,2,.5,.125),10)),6,2)
  ##
  Y_ijk = ds[,c("SampleID_seq",
                "control",
                "dataset_seq",
                "UID_seq",
                "Count_log10")]
  colnames(Y_ijk) <- c("i","j","k","UID","Y")
  controls <- c("A","B","C","D","E","F")
  n_studies = max(Y_ijk$k)
  nk<<-numeric(n_studies)
  for(k in 1 : n_studies)
  {
    n <- nrow(coeffs2[coeffs2$dataset_seq==k,])
    nk[k]<<-n
  }
  n_samples = max(Y_ijk$UID)
  ##Create i/j/k vectors:
  allik_i <<- coeffs2$SampleID_seq
  allik_k <<- coeffs2$dataset_seq
  alljk_j <<- rep(seq(6),n_studies)
  alljk_k <<- rep(seq(n_studies),rep(6,n_studies))
  klist <<- as.list(seq(n_studies))
  Yik_indx <<- 1:length(allik_i)
  Yjk_indx <<- 1:length(alljk_j)

  Yik_list <<- list()

  for (indx in 1:n_samples){
    Yik_list[[indx]] <<- Y_ijk[Y_ijk$i==allik_i[indx] & Y_ijk$k==allik_k[indx],]
  }

  Yjk_list <<- list()

  for (k in 1:n_studies){
    Yjk_list[[k]] <<- list()
    for (j in 1:6){
      Yjk_list[[k]][[j]] <<- Y_ijk[Y_ijk$j==controls[j] & Y_ijk$k==k,]
    }
  }
  ###Randomize Initial Starting values
  ab <- unique(ds[, c("SampleID_seq", "dataset_seq")])
  ab$a <- coeffs2$intercept ##Don't randomize these, since they are sampled first
  ab$b <- coeffs2$slope
  colnames(ab) <- c("i","k","a","b")
  m <- numeric(n_studies*2)
  Sig <- list()
  sigjk <- numeric(n_studies*6)
  sjk <- numeric(n_studies*6)
  for (k in 1:n_studies){
    LB <- min(coeffs2$UID_seq[coeffs2$dataset_seq==k])
    UB <- max(coeffs2$UID_seq[coeffs2$dataset_seq==k])
    a <- coeffs2$intercept[LB:UB]
    b <- coeffs2$slope[LB:UB]
    m[2*k-1] <- mean(a) +runif(1,-3,5)
    m[2*k] <- mean(b) +runif(1,-1,2)
    Sig[[k]]<-matrix(c(var(a),cov(a,b),cov(a,b),var(b)),2,2)*(1/rgamma(1,2,2))
    sigjk[(k*6-5):(k*6)] <- as.vector(tapply(ds$Count_log10[ds$dataset_seq==k],
                                             ds$control[ds$dataset_seq==k],var))*(1/rgamma(1,2,2))
    sjk[(k*6-5):(k*6)] <- as.vector(tapply(ds$residual[ds$dataset_seq==k],
                                           ds$control[ds$dataset_seq==k], mean))*(1/rgamma(1,2,2)) ##runif(6,-.8,.8)
  }
  m <- matrix(m,2,n_studies)
  sigjk <- matrix(sigjk,6,n_studies)
  sjk <- matrix(sjk,6,n_studies)

  mu <- matrix(apply(m,1,mean),2,1) + runif(2,-1,3)
  sig_alpha <- var(m[1,])*(1/rgamma(1,2,2))
  sig_beta <- var(m[2,])*(1/rgamma(1,2,2))
  tj <- matrix(apply(sjk,1,mean),6,1) + runif(6,-.3,.3)
  tau2j <- apply(sjk,1,var)*runif(1,.25,4)

  ###Create DataFrame for Samples
  K <- list()
  K$m <- array(0,c(n_keep,2,n_studies))
  K$mu_alpha <- numeric(n_keep)
  K$mu_beta <- numeric(n_keep)
  K$sig_alpha <- numeric(n_keep)
  K$sig_beta <- numeric(n_keep)
  K$sjk <- array(0,c(n_keep,6,n_studies))
  K$t <- matrix(0,n_keep,6)
  K$tau2j <- matrix(0,n_keep,6)
  K$sigjk <- array(0,c(n_keep,6,n_studies))

  Draws <- data.frame(Index=seq(n_keep),mu_alpha=0,mu_beta=0,
                     sig_alpha=0,sig_beta=0,t1=0,t2=0,t3=0,t4=0,t5=0,t6=0)
  pb = progress::progress_bar$new(total = M)
  for (ii in 1:M){
    pb$tick()
    ab[,3:4]<-t(update_ab(allik_i,allik_k,Yik_indx,m,Sig,sigjk,sjk)) ##.234 secs
    Sig <- lapply(klist,update_Sig,ab1=ab,m=m) ##.016 secs
    K$m[max(c(1, ii-M+n_keep)),,]<-m[,] <- update_m(1:n_studies,ab,sig_alpha,sig_beta,Sig,mu) ##.009 secs
    K$mu_alpha[max(c(1, ii-M+n_keep))]<- mu[1,1]<-update_mu(-1,5,m[1,],sig_alpha, n_studies) ##mu_alpha, .00002 secs
    K$mu_beta[max(c(1, ii-M+n_keep))]<- mu[2,1]<-update_mu(0,4,m[2,],sig_beta, n_studies) ##mu_beta, .00002 secs
    K$sig_alpha[max(c(1, ii-M+n_keep))] <- sig_alpha <- update_sig(.01,m[1,],mu[1,1], n_studies) ### .00002 secs
    K$sig_beta[max(c(1, ii-M+n_keep))] <- sig_beta <- update_sig(.01,m[2,],mu[2,1], n_studies) ### .00002 secs

    sjk[1,]<-update_sjk(1,1:n_studies,tj,tau2j,sigjk,ab,sjk) ### .016 secs
    sjk[2,]<-update_sjk(2,1:n_studies,tj,tau2j,sigjk,ab,sjk) ### .016 secs
    sjk[3,]<-update_sjk(3,1:n_studies,tj,tau2j,sigjk,ab,sjk) ### .016 secs
    sjk[4,]<-update_sjk(4,1:n_studies,tj,tau2j,sigjk,ab,sjk) ### .016 secs
    sjk[5,]<-update_s5k(1:n_studies,sjk) ### .016 secs
    sjk[6,]<- -apply(sjk[1:5,],2,sum)

    K$sjk[max(c(1, ii-M+n_keep)),,] <- sjk

    tj[1,1] <- update_tj(1,-1,1,sjk,tau2j,tj,n_studies) ### .0002 secs
    tj[2,1] <- update_tj(2,-1,1,sjk,tau2j,tj,n_studies) ### .0002 secs
    tj[3,1] <- update_tj(3,-1,1,sjk,tau2j,tj,n_studies) ### .0002 secs
    tj[4,1] <- update_tj(4,-1,1,sjk,tau2j,tj,n_studies) ### .0002 secs
    tj[5,1] <- ((X[6,2]*sum(tj[1:4])-sum(t(X[1:4,2])%*%tj[1:4]))/(X[5,2]-X[6,2])) ### .0002 secs
    tj[6,1] <- (-1)*sum(tj[1:5,1])
    K$t[max(c(1, ii-M+n_keep)), ] = tj
    tau2j[] <- update_tau2j(1:6,.01,sjk,tj,n_studies) ### .0002 secs
    sigjk[,] <- update_sigjk(alljk_j,alljk_k,.01,ab,sjk) ##.017 secs
    Sys.sleep(1/100)
  }
  Draws$mu_alpha = K$mu_alpha
  Draws$mu_beta = K$mu_beta
  Draws$sig_alpha = K$sig_alpha
  Draws$sig_beta = K$sig_beta
  Draws[, c("t1","t2", "t3",
            "t4", "t5", "t6")] = K$t
  return(Draws)
}
