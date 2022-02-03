#skip_on_cran()

source("BW_imp_CH.R")

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

bootstrap_param_est(wide=dataset_NM, b=0, bin=bin, tol = 0, silent = TRUE)

test_that("test for edge cases: bootstrap_param_est part",
          {
            expect_error(bootstrap_param_est(wide=dataset_CH, b=0, bin=bin, tol = 0, silent = TRUE),
                         "wide must be a data.frame.")
            expect_error(bootstrap_param_est(wide=dataset_NM, b=0, bin=as.character(bin), tol = 0, silent = TRUE),
                         "bin must be numeric.")
            expect_error(bootstrap_param_est(wide=dataset_NM, b=0, bin=as.matrix(bin, ncol = 2), tol = 0, silent = TRUE),
                         "bin must be a vector")
            expect_error(bootstrap_param_est(wide=dataset_NM, b=0, bin=rep(1, 30), tol = 0, silent = TRUE))
            expect_error(bootstrap_param_est(wide=dataset_NM, b="0", bin=bin, tol = 0, silent = TRUE),
                         "b must be numeric.")
            expect_error(bootstrap_param_est(wide=dataset_NM, b=c(0,3), bin=bin, tol = 0, silent = TRUE))
            expect_error(bootstrap_param_est(wide=dataset_NM, b=-3, bin=bin, tol = 0, silent = TRUE),
                         "b must be >=0.")
          })


test_that("test for edge cases: impute part",
          {
            expect_error(impute(dataset_CH, m=100, listFormatOut = TRUE, silent = TRUE),
                         "wide must be a data.frame.")
            expect_error(impute(dataset_NM, m=100, listFormatOut = TRUE, silent = TRUE, bin = as.character(bin)),
                         "bin must be numeric.")
            expect_error(impute(dataset_NM, m=100, listFormatOut = TRUE, silent = TRUE, bin = as.matrix(bin, ncol = 2)),
                         "bin must be a vector")
            expect_error(impute(dataset_NM, m=100, listFormatOut = TRUE, silent = TRUE, bin = rep(1, 30)))
            expect_error(impute(dataset_NM, m="100", listFormatOut = TRUE, silent = TRUE),
                         "m must be numeric.")
            expect_error(impute(dataset_NM, m=c(100, 32), listFormatOut = TRUE, silent = TRUE))
            expect_error(impute(dataset_NM, m=1, listFormatOut = TRUE, silent = TRUE),
                         "m must be >=2.")
          })
 
dataset_NM$ref=sample(c(TRUE,FALSE), size=nrow(dataset_NM),replace = TRUE)
test_that("test for edge cases: impute part",
          {
            expect_error(impute_ref(dataset_CH, ref="ref", m=100, listFormatOut = TRUE, silent = TRUE),
                         "wide must be a data.frame.")
            expect_error(impute_ref(dataset_NM, ref="ref", m=100, listFormatOut = TRUE, silent = TRUE, bin = as.character(bin)),
                         "bin must be numeric.")
            expect_error(impute_ref(dataset_NM, ref="ref", m=100, listFormatOut = TRUE, silent = TRUE, bin = as.matrix(bin, ncol = 2)),
                         "bin must be a vector")
            expect_error(impute_ref(dataset_NM, ref="ref", m=100, listFormatOut = TRUE, silent = TRUE, bin = rep(1, 30)))
            expect_error(impute_ref(dataset_NM, ref="ref", m="100", listFormatOut = TRUE, silent = TRUE),
                         "m must be numeric.")
            expect_error(impute_ref(dataset_NM, ref="ref", m=c(100, 32), listFormatOut = TRUE, silent = TRUE))
            expect_error(impute_ref(dataset_NM, ref="ref", m=1, listFormatOut = TRUE, silent = TRUE),
                         "m must be >=2.")
          })
