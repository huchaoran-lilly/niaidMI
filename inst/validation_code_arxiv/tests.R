### Note that NM is programmer 1 and CH is programmer 2

#######################################################################################################

source("../../tests/BW_imp_CH.R")

#######################################################################################################


## bin: a vector with length(dat) - 1 
bin <- c(rep(1, 6), #(Please not data has 28 days)
         rep(2, 7), rep(3, 7), rep(4, 7))

library(niaidMI)
set.seed(2021)
dataset_NM=sim_data(n=200)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]


#Check Model Fit:
fit_CH <- BW_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP)
fit_NM <- bootstrap_param_est(wide=dataset_NM, b=0, bin=bin, tol = 0)$fit


#<<<<--TODO by Chaoran turn this into testthat---->>>>
fit_CH$llk - fit_NM$logLike_new

fit_CH$transition_prob[[1]]-fit_NM$Tran[[1]]
fit_CH$transition_prob[[2]]-fit_NM$Tran[[2]]
fit_CH$transition_prob[[3]]-fit_NM$Tran[[3]]
fit_CH$transition_prob[[4]]-fit_NM$Tran[[4]]

fit_CH$initial_prob-fit_NM$Pri[[1]]

##################################################################
## Check imputation with no stratification
set.seed(2021)
imp_CH <- imputeBS_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, m = 1000)
set.seed(2021)
imp_NM <- impute(dataset_NM, m=1000, listFormatOut = TRUE)


## reformat result to compare results
reform_result_nonstrata <- function(input, input_format = c("NM", "CH")) {
  output <- NA
  if (input_format == "NM") {
    for (i in seq_len(length(input))) {
      output <- rbind(output, as.matrix(input[[i]][, -c(1,2,3)]))
    }
  } else {
    for (i in seq_len(length(input))) {
      output <- rbind(output, as.matrix(input[[i]]))
    }
  }
  
  output[-1, ]
}

imp_NM <- reform_result_nonstrata(imp_NM, "NM")
imp_CH <- reform_result_nonstrata(imp_CH, "CH")


# test transition matrix
estTransMatrx(imp_CH) - estTransMatrx(imp_NM)
# test initial distribution
estInitDist(imp_CH) - estInitDist(imp_NM)




##################################################################
## Check imputation with stratification
set.seed(2021)
imp_strat_CH <- imputeBS_CH(dataset_CH[[1]], dataset_CH[[2]], bin, start_initP, start_tP, m = 100, by = 1)
set.seed(2021)
imp_strat_NM <- impute(dataset_NM, by="strata", m=100, listFormatOut = TRUE)



## reformat result to compare results
reform_result_strata <- function(input, input_format = c("NM", "CH")) {
  output <- NA
  if (input_format == "NM") {
    for (i in seq_len(length(input))) {
      cart <- as.numeric(input[[i]][, 3])
      output <- as.matrix(cbind(cart, input[[i]][, -c(1,2,3)]))
    }
  } else {
    for (i in seq_len(length(input))) {
      output <- as.matrix(input[[i]])
    }
  }
  
  output[-1, ]
}

imp_strat_NM <- reform_result_strata(imp_strat_NM, "NM")
imp_strat_CH <- reform_result_strata(imp_strat_CH, "CH")


# test transition matrix
estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1]) -
estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])

estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1]) -
estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])

# test initial distribution
estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1]) -
estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])

estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1]) -
estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])



