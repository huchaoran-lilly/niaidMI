### Note that NM is programmer 1 and CH is programmer 2

#######################################################################################################

source("BW_imp_CH.R")


## q functions: a n_q_func by 8 matrix with (i, j)th
## element for q(h^obs=i, h=j)
## q_func <- rbind(diag(8), 1)

## bin: a vector with length(dat) - 1 
bin <- c(rep(1, 6), #(Please not data has 28 days)
         rep(2, 7), rep(3, 7), rep(4, 7))

# # ## initial probability of chain
# initP <- c(rep(1, 7) / 7, 0)
# start_initP <- initP
# 
# ## transition probability: a list with length equals
# ## to the number of bins. Each element is a 8*8 matrix.
# tP <- vector('list', length(unique(bin)))
# tP[[1]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
#                  c(rep(0, 7), 1))
# tP[[2]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
#                  c(rep(0, 7), 1))
# tP[[3]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
#                  c(rep(0, 7), 1))
# tP[[4]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
#                  c(rep(0, 7), 1))
# start_tP <- tP






###############################################################
#  Sample Data
###############################################################
library(niaidMI)
dataset_NM=sim_data(n=200)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_initP <- get_start(dataset_CH[[1]][, -1], bin)[[2]]
start_tP <- get_start(dataset_CH[[1]][, -1], bin)[[1]]

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
imp_CH=list()
imp_CH[[1]] <- impute_CH(dataset_CH, q_func, bin, result[[1]], result[[2]]) #todo:Chaoran to replace results with bootstrap results
imp_CH[[2]] <- impute_CH(dataset_CH, q_func, bin, result[[1]], result[[2]]) #todo:Chaoran to replace results with bootstrap results

set.seed(2021)
imp_NM <- impute(dataset_NM, m=2, listFormatOut = TRUE)

#todo: Chaoran to do test that for same


## Check imputation with stratification
set.seed(2021)

#todo: Chaoran to make code for stratified imputation
# Bootstrap + Estimate + impute each strata separately
imp_strat_CH=list()
imp_strat_CH[[1]] <- impute_strat_CH(dataset_CH, q_func, bin, result[[1]], result[[2]])
imp_strat_CH[[2]] <- impute_strat_CH(dataset_CH, q_func, bin, result[[1]], result[[2]])


set.seed(2021)
imp_strat_NM <- impute(dataset_NM, by="strata", m=2, listFormatOut = TRUE)

#todo: Chaoran to make testthat code for comparing stratified imputation
