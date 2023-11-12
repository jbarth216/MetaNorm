#' Estimate mRNA amount for every gene/sample combination
#'
#' This function is similar to RCRnorm with the prior changed
#' to an informative prior infered from the meta analysis. The underlying
#' model and inference procedure is exactly the same as the MetaNorm function.
#' However, it only keeps the most recent draws for each parameter except for
#' global parameters whose posterior means are computed via a rolling mean method.
#' @param dat A list of 4 dataframes
#' @param M Number of MCMC draws
#' @param burn_in Number of burn in MCMC draws
#' @param thin Number of thinning in MCMC draws
#' @param mm The radius
#' @param seed Random seed
#' @return A list containing two items: “GlobalParameters_Draws” (all draws for mu_a, mu_b, sig2_a, sig2_b) and “Posterior_estimates” (list of kappa_hk, kappa_reg posterior estimates)
#' @export
MetaNorm_LowMemory <- function(dat,
                               M=15000,
                               burn_in=5000,
                               thin=1,
                               mm=3,
                               seed=0){

  set.seed(seed)
  tds <- round((M - burn_in)/thin)

  if(burn_in>=M){stop("burn_in must be less than M")}
  if(thin >= (M-burn_in)){stop("thin must be less than the number of post burn-in draws (M-burn_in)")}

  pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125))
  pos_dat = log10(dat$pos_dat + 1)
  n_hk = dim(dat$hk_dat)[1]
  n_reg = dim(dat$reg_dat)[1]
  n_neg = dim(dat$neg_dat)[1]
  n_pos = dim(dat$pos_dat)[1]
  n_patient <- n_samples <- dim(dat$pos_dat)[2]
  X_pos <- matrix(c(rep(1,6),pos_conc),6,2)
  Y_pos <- matrix(log10(unlist(dat$pos_dat)+1),n_pos,n_samples)
  Y_neg <- matrix(log10(unlist(dat$neg_dat)+1),n_neg,n_samples)
  Y_hk <- matrix(log10(unlist(dat$hk_dat)+1),n_hk,n_samples)
  Y_reg <- matrix(log10(unlist(dat$reg_dat)+1),n_reg,n_samples)
  Sig_uns <- solve(outer(1:4,1:4,X_thing,pos_conc)+diag(4))
  D5<-matrix((X_pos[1:4,2]-X_pos[5,2])/(X_pos[5,2]-X_pos[6,2]),4,1)
  D6<-matrix((X_pos[1:4,2]-X_pos[6,2])/(X_pos[5,2]-X_pos[6,2]),4,1)

  Ymat<-rbind(Y_hk,Y_reg)

  ## Create Priors ####
  all_coef = apply(pos_dat, 2, fitWithPosCtrl, pos_conc)

  mu_a_itm = mean(all_coef[1,])
  mu_b_itm = mean(all_coef[2,])

  mu_priors <- c(2.3565,.9518,.0608^2,.0166^2)
  mu_a_mu <- mu_priors[1]
  mu_b_mu <- mu_priors[2]
  sigma2_mu_a <- mu_priors[3]
  sigma2_mu_b <- mu_priors[4]

  hk_RNA = sweep(sweep(log10(dat$hk_dat+1), 2, all_coef[1,], '-'), 2, all_coef[2,], '/')
  reg_RNA = sweep(sweep(log10(dat$reg_dat+1), 2, all_coef[1,], '-'), 2, all_coef[2,], '/')

  ##estimate genes' mean expression level range
  lambda_hk_range = apply(hk_RNA, 1, get_range, mm = mm)
  lambda_reg_range = apply(reg_RNA, 1, get_range, mm = mm)

  #estimate patient effect range by two way ANOVA with patient's regular gene expression level.
  gene = factor(rep(1:n_reg, n_samples))
  patient = factor(rep(1:n_samples, each = n_reg))

  mod = stats::lm(unlist(reg_RNA) ~ patient + gene, contrasts = list(patient = 'contr.sum', gene = 'contr.sum'))

  phi = numeric(n_samples)
  phi[1:(n_samples - 1)] = summary(mod)$coefficients[2:n_samples, 1]
  phi[n_samples] = -sum(phi)

  phi_L = phi - mm * summary(mod)$coefficients[2, 2]
  phi_U = phi + mm * summary(mod)$coefficients[2, 2]


  ##initialize Starting Parameters ####

  sig2_a = stats::runif(1, 0, .01)
  sig2_b = stats::runif(1, 0, .01)
  a = stats::rnorm(n_patient, 2.5, .1)
  mu_a = rnorm(1, 2.5, .1)
  b = stats::rnorm(n_patient, .9, .1)
  mu_b = rnorm(1,.9,.1)
  cc = stats::runif(1, -6, -1)
  phi = stats::rnorm(n_patient, 0, 2)
  phi[n_patient] = -sum(phi[1:(n_patient - 1)])
  kappa_hk = matrix(stats::rnorm(n_hk * n_patient, 0, 1),n_hk,n_patient)
  sig2_kappa_hk = stats::runif(1, 0, 1)
  kappa_reg = matrix(stats::rnorm(n_reg * n_patient, 0, 1),n_reg,n_patient)
  sig2_kappa_reg = stats::runif(1, 0, 1)
  lambda_hk = stats::rnorm(n_hk, 0, 1)
  lambda_reg = stats::rnorm(n_reg, 0, 1)
  d_neg = stats::rnorm(n_neg, 0, .01)
  sig2_dn = stats::runif(1, 0, .1)
  d_pos = stats::rnorm(n_pos, 0, .01)
  sig2_d = stats::runif(1, 0, .1)
  d_hk = stats::rnorm(n_hk, 0, .01)
  d_reg = stats::rnorm(n_reg, 0, .01)
  sig2_n = stats::runif(1, 0, .1)
  sig2_e = stats::runif(1, 0, .1)

  ## Create Output List ####
  Results <- list()

  Results$mu_a <- numeric(M)
  Results$mu_b <- numeric(M)
  Results$sig2_a <- numeric(M)
  Results$sig2_b <- numeric(M)

  keep_a <- matrix(0,1,n_samples)
  keep_b <- matrix(0,1,n_samples)
  keep_cc <- numeric(1)
  keep_d_pos <- matrix(0,1,n_pos)
  keep_d_neg <- matrix(0,1,n_neg)
  keep_d_hk <- matrix(0,1,n_hk)
  keep_d_reg <- matrix(0,1,n_reg)
  keep_phi <- matrix(0,1,n_samples)
  keep_kappa_hk <- matrix(0,n_hk,n_samples)
  keep_kappa_reg <- matrix(0,n_reg,n_samples)
  keep_lambda_hk <- matrix(0,1,n_hk)
  keep_lambda_reg <- matrix(0,1,n_reg)
  keep_sig2_kappa_hk <- numeric(1)
  keep_sig2_kappa_reg <- numeric(1)
  keep_sig2_d <- numeric(1)
  keep_sig2_dn <- numeric(1)
  keep_sig2_e <- numeric(1)
  keep_sig2_n <- numeric(1)

  ##Begin Loop ####
  pb = progress::progress_bar$new(total = M)
  for (i in 1:M){
    pb$tick()
    a <- update_norm_a(Y_pos,X_pos,Y_neg,Y_hk,Y_reg,
                       cc,b,d_pos,d_neg,d_hk,d_reg,mu_a,phi,kappa_hk,kappa_reg,
                       sig2_a,sig2_e,sig2_n,n_pos,n_neg,n_hk,n_reg,n_samples)
    b <- update_norm_b(Y_pos,X_pos,Y_neg,Y_hk,Y_reg,
                       cc,a,d_pos,d_neg,d_hk,d_reg,mu_b,phi,kappa_hk,kappa_reg,
                       sig2_b,sig2_e,sig2_n,n_pos,n_neg,n_hk,n_reg,n_samples)
    cc <- update_norm_cc(Y_neg,a,b,d_neg,sig2_n,L=-6,U=-1,n_neg,n_samples)
    Results$mu_a[i] <- mu_a <- update_norm_mu(a,sig2_a,mu_a_mu,sigma2_mu_a,L=-10,U=10,n_samples)
    Results$mu_b[i] <- mu_b <- update_norm_mu(b,sig2_b,mu_b_mu,sigma2_mu_b,L=-10,U=10,n_samples)

    Y_fixed<-(Ymat - matrix(a,n_hk+n_reg,n_samples,byrow=T))/matrix(b,n_hk+n_reg,n_samples,byrow=T)

    lambda_hk <- apply(Y_fixed[1:n_hk,],1,mean)
    lambda_reg <- apply(Y_fixed[(n_hk+1):(n_hk+n_reg),],1,mean)
    phi <- apply(Y_fixed,2,mean) - mean(Y_fixed)
    kappa_hk <- update_norm_kappa(Y_hk,a,b,phi,d_hk,lambda_hk,
                                  sig2_e,sig2_kappa_hk,n_hk,n_samples)
    kappa_reg <- update_norm_kappa(Y_reg,a,b,phi,d_reg,lambda_reg,
                                   sig2_e,sig2_kappa_reg,n_reg,n_samples)

    for (j in 1:(n_neg-1)){
      d_neg[j] <- update_norm_dneg2(Y_neg,cc,a,b,d_neg,sig2_dn,sig2_n,n_samples,n_neg,j)
    }
    d_neg[n_neg] <- -1*sum(d_neg[-(n_neg)])

    d_pos <- update_norm_dpos3(Y_pos,X_pos,a,b,sig2_d,sig2_e,n_pos,n_samples,Sig_uns,D5,D6)

    d_hk <- update_norm_dhkreg(Y_hk,a,b,phi,kappa_hk,sig2_e,sig2_d,n_hk,n_samples)
    d_reg <- update_norm_dhkreg(Y_reg,a,b,phi,kappa_reg,sig2_e,sig2_d,n_reg,n_samples)
    sig2_kappa_hk <- update_norm_sig2_kappa(kappa_hk,lambda_hk,n_hk,n_samples)
    sig2_kappa_reg <- update_norm_sig2_kappa(kappa_reg,lambda_reg,n_reg,n_samples)
    Results$sig2_a[i] <- sig2_a <- update_norm_sig2_abd(a,mu_a,n_samples)
    Results$sig2_b[i] <- sig2_b <- update_norm_sig2_abd(b,mu_b,n_samples)
    sig2_dn <- update_norm_sig2_abd(d_neg,0,n_neg)
    sig2_d <- update_norm_sig2_abd(d_pos,0,n_pos)
    sig2_n <- update_norm_sig2_n(Y_neg,a,b,cc,d_neg,n_neg,n_samples,eps=.01)
    sig2_e <- update_norm_sig2_e(Y_pos,X_pos,Y_hk,Y_reg,a,b,phi,kappa_hk,kappa_reg,
                                 d_pos,d_hk,d_reg,n_pos,n_hk,n_reg,n_samples,eps=.01)
    if(i>burn_in && i%%thin==0){
      keep_a <- keep_a + a*(1/tds)
      keep_b <- keep_b + b*(1/tds)
      keep_cc <- keep_cc + cc/tds
      keep_lambda_hk <- keep_lambda_hk + lambda_hk/tds
      keep_lambda_reg <- keep_lambda_reg + lambda_reg/tds
      keep_phi <- keep_phi + phi/tds
      keep_kappa_hk <- keep_kappa_hk + kappa_hk/tds
      keep_kappa_reg <- keep_kappa_reg + kappa_reg/tds
      keep_d_neg <- keep_d_neg + d_neg/tds
      keep_d_pos <- keep_d_pos + d_pos/tds
      keep_d_hk <- keep_d_hk + d_hk/tds
      keep_d_reg <- keep_d_reg + d_reg/tds
      keep_sig2_kappa_hk <- keep_sig2_kappa_hk + sig2_kappa_hk/tds
      keep_sig2_kappa_reg <- keep_sig2_kappa_reg + sig2_kappa_reg/tds
      keep_sig2_d <- keep_sig2_d + sig2_d/tds
      keep_sig2_n <- keep_sig2_n + sig2_n/tds
      keep_sig2_e <- keep_sig2_e + sig2_e/tds
      keep_sig2_dn <- keep_sig2_dn + sig2_dn/tds
    }
    Sys.sleep(1/100)
  }

  Output <- list()
  Output$GlobalParameters_Draws <- Results
  Posterior_estimates <- list()
  Posterior_estimates$kappa_hk <- keep_kappa_hk
  Posterior_estimates$kappa_reg <- keep_kappa_reg
  Output$Posterior_estimates <- Posterior_estimates
  return(Output)
}
