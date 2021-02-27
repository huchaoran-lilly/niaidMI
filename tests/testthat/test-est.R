#skip_on_cran()

source("BW_imp_CH.R")

# ## bin: a vector with length(dat) - 1 
bin <- c(rep(1, 6), #(Please not data has 28 days)
         rep(2, 7), rep(3, 7), rep(4, 7))


set.seed(2021)
dataset_NM=sim_data(n=200)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]


#Check Model Fit:
fit_CH <- BW_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, message = FALSE)
fit_NM <- bootstrap_param_est(wide=dataset_NM, b=0, bin=bin, tol = 0, silent = TRUE, maxiter = 100000000000)$fit



#Test Results are equal
test_that("likelihood values are equal between BW and CH",
          {
            expect_equal(fit_CH$llk, fit_NM$logLike_new)
          })

test_that("transition matrix estimate are equal between BW and CH",
          {
            expect_equal(as.vector(fit_CH$transition_prob[[1]]), as.vector(fit_NM$Tran[[1]]))
            expect_equal(as.vector(fit_CH$transition_prob[[2]]), as.vector(fit_NM$Tran[[2]]))
            expect_equal(as.vector(fit_CH$transition_prob[[3]]), as.vector(fit_NM$Tran[[3]]))
            expect_equal(as.vector(fit_CH$transition_prob[[4]]), as.vector(fit_NM$Tran[[4]]))
          })

test_that("initial probability estimate are equal between BW and CH",
          {
            expect_equal(fit_CH$initial_prob, fit_NM$Pri[[1]])
          })

