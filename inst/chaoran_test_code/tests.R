### Note that NM is programmer 1 and CH is programmer 2

#######################################################################################################

source("BW_imp_CH.R")

#######################################################################################################


## bin: a vector with length(dat) - 1 
bin <- c(rep(1, 6), #(Please not data has 28 days)
         rep(2, 7), rep(3, 7), rep(4, 7))

library(niaidMI)
dataset_NM=sim_data(n=100)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]


#Check Model Fit:
fit_CH <- BW_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP)
fit_NM <- bootstrap_param_est(wide=dataset_NM, b=2, bin=bin, tol = 0)[[1]]


#<<<<--TODO by Chaoran turn this into testthat---->>>>
fit_CH$llk - fit_NM$logLike_new

fit_CH$transition_prob[[1]]-fit_NM$Tran[[1]]
fit_CH$transition_prob[[2]]-fit_NM$Tran[[2]]
fit_CH$transition_prob[[3]]-fit_NM$Tran[[3]]
fit_CH$transition_prob[[4]]-fit_NM$Tran[[4]]

fit_CH$initial_prob-fit_NM$Pri[[1]]


## Check imputation with no stratification
set.seed(2021)
imp_CH <- imputeBS_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, m = 2)
set.seed(2021)
imp_NM <- impute(dataset_NM, m=2, listFormatOut = TRUE)

#todo: Chaoran to do test that for same


## Check imputation with stratification
set.seed(2021)
imp_strat_CH <- imputeBS_CH(dataset_CH[[1]], dataset_CH[[2]], bin, start_initP, start_tP, m = 2, by = 1)
set.seed(2021)
imp_strat_NM <- impute(dataset_NM, by="strata", m=2, listFormatOut = TRUE)


## reformat result to compare results
reform_result <- function(input, input_format = c("NM", "CH")) {
  output <- vector('list', length(input))
  if (input_format == "NM") {
    for (i in seq_len(length(input))) {
      output[[i]] <- as.matrix(input[[i]][, -c(1,2,3)])
    }
  } else {
    for (i in seq_len(length(input))) {
      output[[i]] <- as.matrix(input[[i]][, -1])
    }
  }
  
  output
}

imp_strat_NM <- reform_result(imp_strat_NM, "NM")
imp_strat_CH <- reform_result(imp_strat_CH, "CH")

imp_strat_CH[[1]] == imp_strat_NM[[1]]






