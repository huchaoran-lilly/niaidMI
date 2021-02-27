source("BW_imp_CH.R")


bin <- rep(1, 27)

#simulate data
set.seed(2021)
dataset_NM=sim_data(n=200, drop_out_rate = 0.02, sporatic_rate = 0.2)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

#Fit data
start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]
fit_CH <- BW_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, message = FALSE, tol_llk = 1e-10)

initP <- fit_CH[[1]]
tP <- fit_CH[[2]]


#inputs for niaid MI internal function
s=rep(1, nrow(dataset_NM))
Em = get_emission(dataset_NM, paste0("D",1:28))


#inpute both ways
NMimpute <- function(){

  NM <- do.call(rbind, replicate(100,
                           niaidMI:::.impute_one(Pri=list(initP), s=s, 
                                                 Tran=tP, Em=Em, bin=bin,
                                                 fix_const=0),
                           simplify = FALSE),
          )
  NM[-1, ]
}

CHimpute <- function(){
  CH <- do.call(rbind, replicate(100,
                                 impute_CH(dataset_CH[[1]][, -1], 
                                           dataset_CH[[2]], bin = bin, 
                                           initP, tP),
                                 simplify = FALSE),
  )
  CH[-1, ]
}

set.seed(12324)
test1 <- NMimpute()
set.seed(12324)
test2 <- CHimpute()

#test equal
test_that("Imputed values are all equal between BW and CH",
          {
            expect_equal(test1, test2)
          })


