#' Find the simple linear regression coefficients
#'
#' Given the response y and the feature x, this function
#' fits a simple linear regression and returns its coefficients
#' @param y A vector of responses
#' @param x A vector of features
#' @return The coefficients of the fitted simple linear regression
#' @export
fitWithPosCtrl = function(y, x)
{
  mod1 = stats::lm(y ~ x)
  coefs = stats::coef(mod1)
  unname(coefs)
}


#' Get prior range of a uniform distribution
#'
#' This funciton gets a prior range of a uniform distribution
#' This is mainly used for setting up the priors for the Gibbs sampler
#' @param x A vector
#' @param mm Number of standard deviations
#' @return The lower and upper bounds of a uniform distribution
#' @export
get_range = function(x, mm = 5)
{
  c(mean(x)-mm*stats::sd(x), mean(x)+mm*stats::sd(x))
}


#' Get residual from positive control fitted with simple linear regression
#'
#' This function extracts residuals from a series of fitted
#' simple linear regressions
#' @param log_data A matrix of log-transformed data
#' @param RNA_conc A array of RNA concentrations
#' @param coefs A matrix of coefficients of simple linear regressions
#' @return A matrix of residuals
#' @export
get_residual = function(log_dat, RNA_conc, coefs)
{
  log_dat - sweep(sweep(RNA_conc, 2, coefs[2, ], '*'), 2, coefs[1, ], '+')
}


#' Meta analysis data curation
#'
#' This function curates dataset for meta analysis
#' @param dataset This can either be a data.frame or the filename to the dataset
#' @param ... Other paramters to be used by the funtion read.csv
#' @return A dataframe with curated data
#' @export
curate_data = function(dataset, ...)
{
  # We first read in the data
  if(class(dataset) == "character")
  {
    df = read.csv(dataset, ...)
  }else{
    df = dataset
  }
  # The minimum requirement for columns names
  # is that it contains DataSet, RNA, UID, SampleID, and Count
  if(! all(c("DataSet", "RNA", "UID", "SampleID", "Count") %in% colnames(df)))
  {
    stop('DataSet, RNA, UID, SampleID, and Count must be in the column names')
  }
  # There should be only 6 RNA levels and
  # the RNA levels should be 128, 32, 8, 2, 0.5, and 0.125
  rnas = sort(unique(df$RNA), decreasing = TRUE)
  if(length(rnas) != 6)
  {
    stop("There should be 6 unique RNA levels.")
  }
  if(! all(c(128, 32, 8, 2, 0.5, 0.125) == rnas))
  {
    stop('RNA levels should be 128, 32, 8, 2, 0.5, and 0.125')
  }
  # Then we sort the data
  # First by RNA in descending order
  df = df[order(df$RNA, decreasing = TRUE), ]
  # Then by SampleID in ascending order
  df = df[order(df$SampleID, decreasing = FALSE), ]
  # Finally, by DataSet in ascending order
  df = df[order(df$DataSet, decreasing = FALSE), ]
  #################################################
  # Rename DataSet, SampleID so that
  # they are consecutive numbers
  DS_nums = sort(unique(df$DataSet), decreasing = FALSE)
  for(k in 1 : length(DS_nums))
  {
    # We first get the number of sample in the kth study
    n <- length(unique(df$SampleID[df$DataSet==DS_nums[k]]))
    df$SampleID_seq[df$DataSet==DS_nums[k]] <- rep(seq(n),rep(6,n))
    df$dataset_seq[df$DataSet==DS_nums[k]] <- k
  }
  nn <- nrow(df)/6
  df$UID_seq <- rep(seq(nn),rep(6,nn))
  # We add names to RNA levels
  controls <- c("A","B","C","D","E","F")
  df$control = rep(controls, nn)
  # Finally, we add log10 transformed columns
  # to facilitate linear regression
  df$RNA_log10 = log10(df$RNA)
  df$Count_log10 = log10(df$Count+1)
  return(df)
}


#' Find regression coefficients of each patient
#'
#' This function finds regression coefficients of each patient
#' @param df A curated data.frame
#' @return A list that contains the data.frame with residuals and a data.frame containing all regression coefficients
#' @export
find_regression_coefs = function(df)
{
  coeffs2 = unique(df[,c("dataset_seq", "SampleID_seq", "UID_seq")])
  coeffs2$intercept = 0
  coeffs2$slope = 0
  coeffs2$Rsquared = 0
  df$residual = 0
  for(r in 1 : nrow(coeffs2))
  {
    uid = coeffs2$UID_seq[r]
    ind = which(df$UID_seq == uid)
    y = df$Count_log10[ind]
    x = df$RNA_log10[ind]
    reg = summary(lm(y~x))
    df$residual[ind] = reg$residuals
    coeffs2$intercept[r] = reg$coefficients[1,1]
    coeffs2$slope[r] = reg$coefficients[2,1]
    coeffs2$Rsquared[r] = reg$r.squared
  }
  return(list(df=df,
              coeffs2=coeffs2))
}





