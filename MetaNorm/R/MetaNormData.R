#' Datasets of MetaNorm
#'
#' Four .csv files are provided as example data on how to
#' utilize the MetaNorm package.
#'
#' meta_analysis_data.csv contains 13 datasets used for meta-analysis
#' which provided the MetaNorm function with an informative prior.
#'
#' normalization_pos_data.csv, normalization_neg_data.csv,
#' normalization_hk_data.csv, normalization_reg_data.csv contain
#' mRNA measurements on 4 kinds of probes from one study. The data
#' was used by RCRnorm.
#'
#' We seperated the original .RData file into
#' 4 .csv files so that prosepective users can use them as a guidence
#' on how to organize their own data.
#'
#' We also included 3 curated datasets based on meta_analysis_data.csv
#' with which users can compare to make sure that the output for their
#' data is of correct format.
#'
#' The meta_analysis_data_coeff.csv contains our empirical estimates
#' of intercepts and slopes
#'
#' @docType data
#'
"MetaNorm"
