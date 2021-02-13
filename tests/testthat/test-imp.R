skip_on_cran()

source("../BW_imp_CH.R")

## bin: a vector with length(dat) - 1 
# bin <- c(rep(1, 6), #(Please note data has 28 days)
#          rep(2, 7), rep(3, 7), rep(4, 7))
bin <- rep(1, 27)

set.seed(2021)
dataset_NM=sim_data(n=200)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]


##################################################################
## Check imputation with no stratification
set.seed(2021)
imp_CH <- imputeBS_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, m = 100)
set.seed(2021)
imp_NM <- impute(dataset_NM, m=100, listFormatOut = TRUE, silent = TRUE)


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

tol_nonstrata <- 0.01

test_that("non-strata imputation functions works the same between NM and CH",
          {
            # test transition matrix
            expect_true(all((abs(estTransMatrx(imp_CH) - estTransMatrx(imp_NM))) <= tol_nonstrata))
            # test initial distribution
            expect_true(all((abs(estInitDist(imp_CH) - estInitDist(imp_NM))) <= tol_nonstrata))
          })



##################################################################
## Check imputation with stratification
set.seed(2021)
imp_strat_CH <- imputeBS_CH(dataset_CH[[1]], dataset_CH[[2]], bin, start_initP, start_tP, m = 100, by = 1)
set.seed(2021)
imp_strat_NM <- impute(dataset_NM, by="strata", m=100, listFormatOut = TRUE, silent = TRUE)



## reformat result to compare results
reform_result_strata <- function(input, input_format = c("NM", "CH")) {
  output <- NA
  if (input_format == "NM") {
    for (i in seq_len(length(input))) {
      cart <- as.numeric(as.factor(input[[i]][, 3]))
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

tol_strata <- 0.05

# test transition matrix
trn1 <- all(abs(estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1]) - estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])) <= tol_strata)

trn2 <- all(abs(estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1]) - estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])) <= tol_strata)

# test initial distribution
int1 <- all(abs(estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1]) - estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])) <= tol_strata)

int2 <- all(abs(estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1]) - estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])) <= tol_strata)


test_that("strata imputation functions work the same between NM and CH",
          {
            expect_true(trn1)
            expect_true(trn2)
            expect_true(int1)
            expect_true(int2)
          })














