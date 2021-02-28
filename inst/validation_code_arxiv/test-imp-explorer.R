library(niaidMI)
source("BW_imp_CH_debug.R")

## bin: a vector with length(dat) - 1 
# bin <- c(rep(1, 6), #(Please note data has 28 days)
#          rep(2, 7), rep(3, 7), rep(4, 7))
bin <- rep(1, 27)

set.seed(2021)
#dataset_NM=sim_data(n=200)
dataset_NM=sim_data(n=200, drop_out_rate = 0.02, sporatic_rate = 0.2)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]

#setSessionTimeLimit(cpu = Inf, elapsed = Inf)
##################################################################
## Check imputation with no stratification
set.seed(2021)
imp_CH <- imputeBS_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, m = 50,
                      tol_llk = 1e-10, dataNM = dataset_NM)

set.seed(2021)
imp_NM <- impute(dataset_NM, m=50, listFormatOut = TRUE, silent = TRUE, tol = 1e-10, maxiter = 100000)

# set.seed(2021)
# imp_NM <- impute_debug(dataset_NM, m=50, listFormatOut = TRUE, silent = TRUE, tol = 1e-10, maxiter = 100000)

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

imp_NMt <- reform_result_nonstrata(imp_NM, "NM")
imp_CHt <- reform_result_nonstrata(imp_CH, "CH")

table(imp_NMt == imp_CHt)


# tol_nonstrata <- 0.0001
# 
# test_that("non-strata imputation functions works the same between NM and CH",
#           {
#             # test transition matrix
#             expect_true(all((abs(estTransMatrx(imp_CH) - estTransMatrx(imp_NM))) <= tol_nonstrata))
#             # test initial distribution
#             expect_true(all((abs(estInitDist(imp_CH) - estInitDist(imp_NM))) <= tol_nonstrata))
#           })



##################################################################
## Check imputation with stratification
set.seed(2021)
imp_strat_CH <- imputeBS_CH(dataset_CH[[1]], dataset_CH[[2]], bin, start_initP, start_tP,
                            m = 50, by = 1,tol_llk = 1e-10, dataNM = dataset_NM)
set.seed(2021)
imp_strat_NM <- impute(dataset_NM, by="strata", m=50, listFormatOut = TRUE, silent = TRUE,
                       tol = 1e-10, maxiter = 100000)



## reformat result to compare results
reform_result_strata <- function(input, input_format = c("NM", "CH")) {
  output <- NA
  if (input_format == "NM") {
    for (i in seq_len(length(input))) {
      cart <- as.numeric(as.factor(input[[i]][, 3]))
      output <- rbind(output, as.matrix(cbind(cart, input[[i]][, -c(1,2,3)])))
    }
  } else {
    for (i in seq_len(length(input))) {
      output <- rbind(output, as.matrix(input[[i]]))
    }
  }
  
  output[-1, ]
}

imp_strat_NMt <- reform_result_strata(imp_strat_NM, "NM")
imp_strat_CHt <- reform_result_strata(imp_strat_CH, "CH")

table(imp_strat_NMt == imp_strat_CHt)

table(dataset_CH[[1]][, -1] %in% 1:8)

(191+6) / (1825 * 50 * 2)

# tol_strata <- 0.05
# 
# # test transition matrix
# trn1 <- all(abs(estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1]) - estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])) <= tol_strata)
# 
# trn2 <- all(abs(estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1]) - estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])) <= tol_strata)
# 
# # test initial distribution
# int1 <- all(abs(estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1]) - estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])) <= tol_strata)
# 
# int2 <- all(abs(estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1]) - estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])) <= tol_strata)
# 
# 
# test_that("strata imputation functions work the same between NM and CH",
#           {
#             expect_true(trn1)
#             expect_true(trn2)
#             expect_true(int1)
#             expect_true(int2)
#           })














