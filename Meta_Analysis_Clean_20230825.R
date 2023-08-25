###This file contains code to run the meta analysis MCMC
###Also includes diagnostics and plot checks

##NanoString prior Meta-Analysis
  library(MASS)
  #install.packages("rbenchmark")
  library(rbenchmark)
  library(Rcpp)
  library(truncnorm)
  #setwd("C:\\Users\\jbart\\Documents\\Summer 2021 RA work\\RCRNorm")
  sourceCpp("armadillo.cpp")
  ##write.csv(coeffs2,"coeffs2.csv")
All_clean2 <- read.csv("All_clean2.csv")
coeffs2 <- read.csv("coeffs2.csv")
  ds <- All_clean2 ##See RCRnorm_exp.R. coeffs2 also comes from this script
##We are sequentializing i and k here (since some datasets and samples are removed)

  nk <- numeric(13)
  for (k in 1:13){
    n <- nrow(coeffs2[coeffs2$ds==DS_nums[k],])
    ds$SampleID_seq[ds$DataSet==DS_nums[k]] <- rep(seq(n),rep(6,n))
    ds$dataset_seq[ds$DataSet==DS_nums[k]] <- k
    coeffs2$SampleID_seq[coeffs2$ds==DS_nums[k]] <- seq(n)
    coeffs2$ds_seq[coeffs2$ds==DS_nums[k]] <- k
    nk[k] <- n
  }
  nn <- nrow(ds)/6
  ds$UID_seq <- rep(seq(nn),rep(6,nn))
  coeffs2$UID_seq <- seq(nn)


Y_ijk <- ds[,c(10,2,11,12,7)]
colnames(Y_ijk) <- c("i","j","k","UID","Y")
controls <- c("A","B","C","D","E","F")
Y_list <- list()

##Create i/j/k vectors:
allik_i <- coeffs2$SampleID_seq
allik_k <- coeffs2$ds_seq
alljk_j <- rep(seq(6),13)
alljk_k <- rep(seq(13),rep(6,13))
klist <- as.list(seq(13))
Yik_indx <- 1:length(allik_i)
Yjk_indx <- 1:length(alljk_j)

Yik_list <- list()

for (indx in 1:9929){
  Yik_list[[indx]] <-Y_ijk[Y_ijk$i==allik_i[indx] & Y_ijk$k==allik_k[indx],]
}

Yjk_list <- list()

for (k in 1:13){
  Yjk_list[[k]] <-list()
  for (j in 1:6){
   Yjk_list[[k]][[j]] <- Y_ijk[Y_ijk$j==controls[j] & Y_ijk$k==k,]
  }
}

##Some items that won't change##
I2 <- matrix(c(1,0,0,1),2,2)
X <- matrix(c(rep(1,6),log(c(128,32,8,2,.5,.125),10)),6,2)

##Conditionals as functions

update_ab<-Vectorize(function(i,k,indx,m,Sig,sigjk,sjk){
  m <- m[,k]
  Sig <- Sig[[k]]
  sigjk <- 1/sigjk[,k]
  Y <- Yik_list[[indx]]$Y
  s <- sjk[,k]
  Sig1 <- Sig1_cpp(Sig,X,diag(sigjk))
  mu <- mu_cpp(Sig1,Sig,X,diag(sigjk),m,Y,s)
  rMvNorm(1,mu,Sig1)
},vectorize.args=c("i","k","indx"))

# update_ab2 <- function(m,Sig,sigj,Y,s){
#   mvrnorm(1,solve(solve(Sig) + t(X)%*%diag(sigj)%*%X)%*%(solve(Sig)%*%m+t(X)%*%diag(sigj)%*%(Y-s))
#           ,solve(solve(Sig) + t(X)%*%diag(sigj)%*%X))
# }
#
# benchmark(update_ab(c(1,1),matrix(c(1,0,0,1),2,2),rep(1,6),c(1,2,3,4,5,6),c(1,-1,2,2,1,1)),
#           update_ab2(c(1,1),matrix(c(1,0,0,1),2,2),rep(1,6),c(1,2,3,4,5,6),c(1,-1,2,2,1,1)))
##Clearly 1 is better

##Sampling from a Wish with inv. of Rate Matrix, then inverting
update_Sig <- function(k,ab1,m){
  mu_rep<-matrix(rep(m[,k],nk[k]),nk[k],2,byrow=T)
  Phi<-as.matrix(t(ab1[ab1$k==k,3:4]-mu_rep))%*%as.matrix((ab1[ab1$k==k,3:4]-mu_rep))
  RateMat<-inverse(matrix(I2+Phi,2,2))
  inverse(matrix(rWishart(1,3+nk[k],RateMat),2,2))
}


update_m <- Vectorize(function(k,ab,sig_alpha,sig_beta,Sig,mu){
  theta <- t(as.matrix(ab[ab$k==k,3:4]))
  Sig_m <- matrix(c(sig_alpha,0,0,sig_beta),2,2)
  Sig_m_inv <- solve(Sig_m)
  Sig_inv <- solve(Sig[[k]])
  Var <- solve(Sig_m_inv+nk[k]*Sig_inv)
  Mean <- Var%*%(Sig_m_inv%*%mu+apply(Sig_inv%*%theta,1,sum))
  mvrnorm(1,Mean,Var)
},vectorize.args="k")

###
#sig_alpha <- 1
#sig_beta <- 1
#Sig <- I2
#n <- 4
#mu <- matrix(c(2.5,.1),2,1)
#benchmark(update_m(a,b,sig_alpha,sig_beta,Sig,n,mu),
#          update_m1(a,b,sig_alpha,sig_beta,Sig,n,mu))
###

update_mu <- function(L,U,alpha,sig_alpha,K=13){
  if (L >= U){stop("Lower bound cannot be larger than Upper bound")}
  Mean <- (1/K)*sum(alpha)
  Sd <- ((1/K)*sig_alpha)^.5
  rtruncnorm(1,L,U,Mean,Sd)
}

#update_mu_alpha(-.2,.2,1,1,1)
#Note that this can be used to update both alpha and beta

update_sig <- function(eps,alpha,mu,K=13){
  a <- eps + (K/2)
  b <- eps + .5*sum((alpha-mu)^2)
  1/rgamma(1,a,b)
}

#Note that this can be used to update both sig_alpha and sig_beta

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

update_s5k <- Vectorize(function(k,sjk){
  s <- sum(sjk[1:4,k])
  Xs <- sum(t(X[1:4,2])%*%sjk[1:4,k])
  (X[6,2]*s - Xs)/(X[5,2]-X[6,2])
},vectorize.args="k")

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

update_tau2j <- Vectorize(function(j,eps,sjk,tj,K=13){
  a <- K/2 + eps
  b <- eps + .5*sum((sjk[j,]-tj[j])^2)
  1/rgamma(1,a,b)
},vectorize.args="j")

update_sigjk <- Vectorize(function(j,k,eps,ab,sjk){
  Y <- Yjk_list[[k]][[j]]$Y
  a <- eps + nk[k]/2
  b <- eps + .5*sum((Y-X[j,]%*%t(as.matrix(ab[ab$k==k,3:4]))-sjk[j,k])^2)
  1/rgamma(1,a,b)
},vectorize.args=c("j","k"))



##Loop
Run_Chain <- function(M){
  ###Randomize Initial Starting values
  ab <- unique(ds[c(10,11)])
  ab$a <- coeffs2$intercept ##Don't randomize these, since they are sampled first
  ab$b <- coeffs2$slope 
  colnames(ab) <- c("i","k","a","b")
  m <- numeric(13*2)
  Sig <- list()
  sigjk <- numeric(13*6)
  sjk <- numeric(13*6)
  for (k in 1:13){
    LB <- min(coeffs2$UID_seq[coeffs2$ds==DS_nums[k]])
    UB <- max(coeffs2$UID_seq[coeffs2$ds==DS_nums[k]])
    a <- coeffs2$intercept[LB:UB]
    b <- coeffs2$slope[LB:UB]
    m[2*k-1] <- mean(a) +runif(1,-3,5)
    m[2*k] <- mean(b) +runif(1,-1,2)
    Sig[[k]]<-matrix(c(var(a),cov(a,b),cov(a,b),var(b)),2,2)*(1/rgamma(1,2,2))
    sigjk[(k*6-5):(k*6)] <- as.vector(tapply(ds$Count_log10[ds$DataSet==DS_nums[k]],
                                             ds$control[ds$DataSet==DS_nums[k]],var))*(1/rgamma(1,2,2))
    sjk[(k*6-5):(k*6)] <- as.vector(tapply(ds$residual[ds$DataSet==DS_nums[k]],
                                           ds$control[ds$DataSet==DS_nums[k]], mean))*(1/rgamma(1,2,2)) ##runif(6,-.8,.8)
  }
  m <- matrix(m,2,13)
  sigjk <- matrix(sigjk,6,13)
  sjk <- matrix(sjk,6,13)
  
  mu <- matrix(apply(m,1,mean),2,1) + runif(2,-1,3)
  sig_alpha <- var(m[1,])*(1/rgamma(1,2,2))
  sig_beta <- var(m[2,])*(1/rgamma(1,2,2))
  tj <- matrix(apply(sjk,1,mean),6,1) + runif(6,-.3,.3)
  tau2j <- apply(sjk,1,var)*runif(1,.25,4)  

###Create DataFrame for Samples
K <- list()
K$m <- array(0,c(M,2,13))
K$mu_alpha <- numeric(M)
K$mu_beta <- numeric(M)
K$sig_alpha <- numeric(M)
K$sig_beta <- numeric(M)
K$sjk <- array(0,c(M,6,13))
K$t <- matrix(0,M,6)
K$tau2j <- matrix(0,M,6)
K$sigjk <- array(0,c(M,6,13))

#Draws <- data.frame(Index=seq(M),mu_alpha=0,mu_beta=0,
#                    sig_alpha=0,sig_beta=0,t1=0,t2=0,t3=0,t4=0,t5=0,t6=0)
for (ii in 1:M){
  ab[,3:4]<-t(update_ab(allik_i,allik_k,Yik_indx,m,Sig,sigjk,sjk)) ##.234 secs
  Sig <- lapply(klist,update_Sig,ab1=ab,m=m) ##.016 secs
  K$m[ii,,]<-m[,] <- update_m(1:13,ab,sig_alpha,sig_beta,Sig,mu) ##.009 secs
  K$mu_alpha[ii]<- mu[1,1]<-update_mu(-1,5,m[1,],sig_alpha) ##mu_alpha, .00002 secs
  K$mu_beta[ii]<- mu[2,1]<-update_mu(0,4,m[2,],sig_beta) ##mu_beta, .00002 secs
  K$sig_alpha[ii] <- sig_alpha <- update_sig(.01,m[1,],mu[1,1]) ### .00002 secs
  K$sig_beta[ii] <- sig_beta <- update_sig(.01,m[2,],mu[2,1]) ### .00002 secs
  
  sjk[1,]<-update_sjk(1,1:13,tj,tau2j,sigjk,ab,sjk) ### .016 secs
  sjk[2,]<-update_sjk(2,1:13,tj,tau2j,sigjk,ab,sjk) ### .016 secs
  sjk[3,]<-update_sjk(3,1:13,tj,tau2j,sigjk,ab,sjk) ### .016 secs
  sjk[4,]<-update_sjk(4,1:13,tj,tau2j,sigjk,ab,sjk) ### .016 secs
  sjk[5,]<-update_s5k(1:13,sjk) ### .016 secs
  sjk[6,]<- -apply(sjk[1:5,],2,sum)
  
  K$sjk[ii,,] <- sjk
  
   tj[1,1] <- update_tj(1,-1,1,sjk,tau2j,tj) ### .0002 secs
   tj[2,1] <- update_tj(2,-1,1,sjk,tau2j,tj) ### .0002 secs
   tj[3,1] <- update_tj(3,-1,1,sjk,tau2j,tj) ### .0002 secs
   tj[4,1] <- update_tj(4,-1,1,sjk,tau2j,tj) ### .0002 secs
   tj[5,1] <- ((X[6,2]*sum(tj[1:4])-sum(t(X[1:4,2])%*%tj[1:4]))/(X[5,2]-X[6,2])) ### .0002 secs
   tj[6,1] <- (-1)*sum(tj[1:5,1])
  
  tau2j[] <- update_tau2j(1:6,.01,sjk,tj) ### .0002 secs
  sigjk[,] <- update_sigjk(alljk_j,alljk_k,.01,ab,sjk) ##.017 secs
}
Draws
}

##Add beeps
library(beepr)

ptm1 <- proc.time()
Res1 <- Run_Chain(12000)

ptm2 <- proc.time()
ptm2-ptm1
beep(1)

Res2 <- Run_Chain(12000)

ptm3 <- proc.time()
ptm3-ptm2
beep(1)

Res3 <- Run_Chain(12000)

ptm4 <- proc.time()
ptm4-ptm3
beep(1)

Res4 <- Run_Chain(12000)

proc.time()-ptm4
proc.time()-ptm1
beep(2)

##Diagnostics

##Res1_keep <- Res1
##Res2_keep <- Res2
##Res3_keep <- Res3
##Res4_keep <- Res4

####Load Results####
# data_read<-read.csv("MCMC_results_20211009.csv")
# Res1 <- data_read[data_read$chain==1,2:6]
# Res2 <- data_read[data_read$chain==2,2:6]
# Res3 <- data_read[data_read$chain==3,2:6]
# Res4 <- data_read[data_read$chain==4,2:6]
####End

x11()
par(mfrow=c(2,2))
plot(Res1$Index,Res1$mu_alpha,type="l",main="mu_alpha",col="red")
plot(Res1$Index,Res1$mu_beta,type="l",col="blue")
plot(Res1$Index,Res1$sig_alpha,type="l",col="green")
plot(Res1$Index,Res1$sig_beta,type="l",col="orange")

mu_alph_min <- min(Res1$mu_alpha,Res2$mu_alpha,Res3$mu_alpha,Res4$mu_alpha)
mu_alph_max <- max(Res1$mu_alpha,Res2$mu_alpha,Res3$mu_alpha,Res4$mu_alpha)


x11()
par(mfrow=c(2,2))
plot(Res1$Index,Res1$mu_alpha,type="l",main="mu_alpha - chain 1",col="red",ylim=c(mu_alph_min,mu_alph_max))
plot(Res2$Index,Res2$mu_alpha,type="l",main="mu_alpha - chain 2",col="blue",ylim=c(mu_alph_min,mu_alph_max))
plot(Res3$Index,Res3$mu_alpha,type="l",main="mu_alpha - chain 3",col="green",ylim=c(mu_alph_min,mu_alph_max))
plot(Res4$Index,Res4$mu_alpha,type="l",main="mu_alpha - chain 4",col="orange",ylim=c(mu_alph_min,mu_alph_max))

mu_bet_min <- min(Res1$mu_beta,Res2$mu_beta,Res3$mu_beta,Res4$mu_beta)
mu_bet_max <- max(Res1$mu_beta,Res2$mu_beta,Res3$mu_beta,Res4$mu_beta)


x11()
par(mfrow=c(2,2))
plot(Res1$Index,Res1$mu_beta,type="l",main="mu_beta - chain 1",col="red",ylim=c(mu_bet_min,mu_bet_max))
plot(Res2$Index,Res2$mu_beta,type="l",main="mu_beta - chain 2",col="blue",ylim=c(mu_bet_min,mu_bet_max))
plot(Res3$Index,Res3$mu_beta,type="l",main="mu_beta - chain 3",col="green",ylim=c(mu_bet_min,mu_bet_max))
plot(Res4$Index,Res4$mu_beta,type="l",main="mu_beta - chain 4",col="orange",ylim=c(mu_bet_min,mu_bet_max))

sig_alph_min <- min(Res1$sig_alpha,Res2$sig_alpha,Res3$sig_alpha,Res4$sig_alpha)
sig_alph_max <- max(Res1$sig_alpha,Res2$sig_alpha,Res3$sig_alpha,Res4$sig_alpha)


x11()
par(mfrow=c(2,2))
plot(Res1$Index,Res1$sig_alpha,type="l",main="sig_alpha - chain 1",col="red",ylim=c(sig_alph_min,sig_alph_max))
plot(Res2$Index,Res2$sig_alpha,type="l",main="sig_alpha - chain 2",col="blue",ylim=c(sig_alph_min,sig_alph_max))
plot(Res3$Index,Res3$sig_alpha,type="l",main="sig_alpha - chain 3",col="green",ylim=c(sig_alph_min,sig_alph_max))
plot(Res4$Index,Res4$sig_alpha,type="l",main="sig_alpha - chain 4",col="orange",ylim=c(sig_alph_min,sig_alph_max))

sig_beta_min <- min(Res1$sig_beta,Res2$sig_beta,Res3$sig_beta,Res4$sig_beta)
sig_beta_max <- max(Res1$sig_beta,Res2$sig_beta,Res3$sig_beta,Res4$sig_beta)

x11()
par(mfrow=c(2,2))
plot(Res1$Index,Res1$sig_beta,type="l",main="sig_beta - chain 1",col="red", ylim=c(sig_beta_min,sig_beta_max))
plot(Res2$Index,Res2$sig_beta,type="l",main="sig_beta - chain 2",col="blue",ylim=c(sig_beta_min,sig_beta_max))
plot(Res3$Index,Res3$sig_beta,type="l",main="sig_beta - chain 3",col="green",ylim=c(sig_beta_min,sig_beta_max))
plot(Res4$Index,Res4$sig_beta,type="l",main="sig_beta - chain 4",col="orange",ylim=c(sig_beta_min,sig_beta_max))

##t trace plots###
x11()
par(mfrow=c(2,3))
plot(Res1$Index,Res1$t1,type="l",main="t1")
plot(Res1$Index,Res1$t2,type="l",main="t2")
plot(Res1$Index,Res1$t3,type="l",main="t3")
plot(Res1$Index,Res1$t4,type="l",main="t4")
plot(Res1$Index,Res1$t5,type="l",main="t5")
plot(Res1$Index,Res1$t6,type="l",main="t6")

corvec1 <- numeric(300)

  for (i in 1:300){
    corvec1[i]<-cor(Res1$t1[3001:11700],
                    Res1$t1[(3001+i):(11700+i)])
  }


x11()
plot(1:300,corvec1,xlab="lag",ylab="autocorrelation", 
     main="t1 autocorrelation post burn-in (chain1)",type="l",ylim=c(-.1,1))
abline(h=0,col="red",lty=2)

corvec1 <- numeric(300)

  for (i in 1:300){
    corvec1[i]<-cor(samples$d_pos[5001:7700,1],
                    samples$d_pos[(5001+i):(7700+i),1])
  }

x11()
plot(1:300,corvec1,xlab="lag",ylab="autocorrelation", 
     main="d_pos autocorrelation post burn-in (chain1)",type="l",ylim=c(-.1,1))
abline(h=0,col="red",lty=2)


##Correlation
##mu_alpha

corvec1 <- numeric(50)
corvec2 <- numeric(50)
corvec3 <- numeric(50)
corvec4 <- numeric(50)

for (i in 1:50){
  for (i in 1:50){
    corvec1[i]<-cor(Res1$mu_alpha[3001:11950],
                    Res1$mu_alpha[(3001+i):(11950+i)])
    corvec2[i]<-cor(Res2$mu_alpha[3001:11950],
                    Res2$mu_alpha[(3001+i):(11950+i)])
    corvec3[i]<-cor(Res3$mu_alpha[3001:11950],
                    Res3$mu_alpha[(3001+i):(11950+i)])
    corvec4[i]<-cor(Res4$mu_alpha[3001:11950],
                    Res4$mu_alpha[(3001+i):(11950+i)])
  }
}


x11()
par(mfrow=c(2,2))
plot(1:50,corvec1,xlab="lag",ylab="autocorrelation", 
     main="mu_alpha autocorrelation post burn-in (chain1)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec2,xlab="lag",ylab="autocorrelation", 
     main="mu_alpha autocorrelation post burn-in (chain2)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec3,xlab="lag",ylab="autocorrelation", 
     main="mu_alpha autocorrelation post burn-in (chain3)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec4,xlab="lag",ylab="autocorrelation", 
     main="mu_alpha autocorrelation post burn-in (chain4)",type="b")
abline(h=0,col="red",lty=2)
##mu_beta

corvec1 <- numeric(50)
corvec2 <- numeric(50)
corvec3 <- numeric(50)
corvec4 <- numeric(50)

for (i in 1:50){
  for (i in 1:50){
    corvec1[i]<-cor(Res1$mu_beta[3001:11950],
                    Res1$mu_beta[(3001+i):(11950+i)])
    corvec2[i]<-cor(Res2$mu_beta[3001:11950],
                    Res2$mu_beta[(3001+i):(11950+i)])
    corvec3[i]<-cor(Res3$mu_beta[3001:11950],
                    Res3$mu_beta[(3001+i):(11950+i)])
    corvec4[i]<-cor(Res4$mu_beta[3001:11950],
                    Res4$mu_beta[(3001+i):(11950+i)])
  }
}

x11()
par(mfrow=c(2,2))
plot(1:50,corvec1,xlab="lag",ylab="autocorrelation", 
     main="mu_beta autocorrelation post burn-in (chain1)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec2,xlab="lag",ylab="autocorrelation", 
     main="mu_beta autocorrelation post burn-in (chain2)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec3,xlab="lag",ylab="autocorrelation", 
     main="mu_beta autocorrelation post burn-in (chain3)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec4,xlab="lag",ylab="autocorrelation", 
     main="mu_beta autocorrelation post burn-in (chain4)",type="b")
abline(h=0,col="red",lty=2)
##sig_alpha

corvec1 <- numeric(50)
corvec2 <- numeric(50)
corvec3 <- numeric(50)
corvec4 <- numeric(50)

for (i in 1:50){
  corvec1[i]<-cor(Res1$sig_alpha[3001:11950],
                  Res1$sig_alpha[(3001+i):(11950+i)])
  corvec2[i]<-cor(Res2$sig_alpha[3001:11950],
                  Res2$sig_alpha[(3001+i):(11950+i)])
  corvec3[i]<-cor(Res3$sig_alpha[3001:11950],
                  Res3$sig_alpha[(3001+i):(11950+i)])
  corvec4[i]<-cor(Res4$sig_alpha[3001:11950],
                  Res4$sig_alpha[(3001+i):(11950+i)])
}

x11()
par(mfrow=c(2,2))
plot(1:50,corvec1,xlab="lag",ylab="autocorrelation", 
     main="sig_alpha autocorrelation post burn-in (chain1)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec2,xlab="lag",ylab="autocorrelation", 
     main="sig_alpha autocorrelation post burn-in (chain2)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec3,xlab="lag",ylab="autocorrelation", 
     main="sig_alpha autocorrelation post burn-in (chain3)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec4,xlab="lag",ylab="autocorrelation", 
     main="sig_alpha autocorrelation post burn-in (chain4)",type="b")
abline(h=0,col="red",lty=2)

##sig_beta

corvec1 <- numeric(50)
corvec2 <- numeric(50)
corvec3 <- numeric(50)
corvec4 <- numeric(50)

for (i in 1:50){
  corvec1[i]<-cor(Res1$sig_beta[3001:11950],
                  Res1$sig_beta[(3001+i):(11950+i)])
  corvec2[i]<-cor(Res2$sig_beta[3001:11950],
                  Res2$sig_beta[(3001+i):(11950+i)])
  corvec3[i]<-cor(Res3$sig_beta[3001:11950],
                  Res3$sig_beta[(3001+i):(11950+i)])
  corvec4[i]<-cor(Res4$sig_beta[3001:11950],
                  Res4$sig_beta[(3001+i):(11950+i)])
}

x11()
par(mfrow=c(2,2))
plot(1:50,corvec1,xlab="lag",ylab="autocorrelation", 
     main="sig_beta autocorrelation post burn-in (chain1)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec2,xlab="lag",ylab="autocorrelation", 
     main="sig_beta autocorrelation post burn-in (chain2)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec3,xlab="lag",ylab="autocorrelation", 
     main="sig_beta autocorrelation post burn-in (chain3)",type="b")
abline(h=0,col="red",lty=2)
plot(1:50,corvec4,xlab="lag",ylab="autocorrelation", 
     main="sig_beta autocorrelation post burn-in (chain4)",type="b")
abline(h=0,col="red",lty=2)

##Chain Statistics post burn-in


chain1 <- Res1[3001:12000,]
chain2 <- Res2[3001:12000,]
chain3 <- Res3[3001:12000,]
chain4 <- Res4[3001:12000,]


chain1$chain <- 1
chain2$chain <- 2
chain3$chain <- 3
chain4$chain <- 4

Res1$chain <- 1
Res2$chain <- 2
Res3$chain <- 3
Res4$chain <- 4

chains_all <- rbind(Res1,Res2,Res3,Res4)
chains_all$burn_in <- 0
chains_all$burn_in[chains_all$Index<=3000] <- 1

chains<-rbind(chain1,chain2,chain3,chain4)
x11()
par(mfrow=c(2,2))
boxplot(chains$mu_alpha~chains$chain,main="mu_alpha",xlab="chain",ylab="mu_alpha")
boxplot(chains$mu_beta~chains$chain,main="mu_beta",xlab="chain",ylab="mu_beta")
boxplot(chains$sig_alpha~chains$chain,main="sig_alpha",xlab="chain",ylab="sig_alpha")
boxplot(chains$sig_beta~chains$chain,main="sig_beta",xlab="chain",ylab="sig_beta")

x11()
par(mfrow=c(2,2))
hist(chains$mu_alpha,main="mu_alpha density")
hist(chains$mu_beta,main="mu_beta density")
hist(chains$sig_alpha,main="sig_alpha density")
hist(chains$sig_beta,main="sig_beta density")

#write.csv(chains_all,"MCMC_results_20211009.csv")


#####################################
## Diagnostics and Post. Estimates ##
#####################################

##Thin Chains (every other)
chains_thinned <- chains[chains$Index%%2==1,]

##1) 5-num summary

five_number_summary<-t(apply(chains_thinned[,2:5],2,fivenum))
means <- t(t(apply(chains_thinned[,2:5],2,mean)))
SDs <- t(t(apply(chains_thinned[,2:5],2,sd)))
table1<-cbind(five_number_summary,means,SDs)
c("Minimum","Quarter 1","Median","Quarter 3", 
  "Maximum", "Mean", "Standard Deviation") -> colnames(table1)
write.csv(table1,"Table1.csv")

##2) Gelman Rubin Diagnostics
library(coda)
draws1 <- mcmc(Res1[,2:5])
draws2 <- mcmc(Res2[,2:5])
draws3 <- mcmc(Res3[,2:5])
draws4 <- mcmc(Res4[,2:5])
drawlist <- mcmc.list(list(draws1,draws2,draws3,draws4))
##drawlist <- mcmc.list(list(draws1,draws2))
gelman.diag(drawlist)
x11()
gelman.plot(drawlist)

