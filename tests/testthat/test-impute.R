#skip_on_cran()

source("BW_imp_CH.R")

.impute_one_debug <-
  function(Pri, s, Tran, Em, bin, n_bin=max(bin),
           n=dim(Em)[1], n_day=dim(Em)[2], n_states=8, n_strata=max(s), fix_const=1E-6) {
    
    # 
    # Pri=lapply(Pri,function(x) {
    #   x=x+fix_const
    #   x/sum(x)})
    # Tran=lapply(Tran,function(x) {
    #   x=x+fix_const
    #   x/rowSums(x)})
    
    fb <- niaidMI:::.forward_backward(Pri=Pri, s=s, Tran=Tran, Em=Em, bin=bin,
                                      n_day=n_day, n_states=n_states, n_bin=n_bin, n_strata=n_strata)
    
    ret=matrix(0, n, n_day)
    
    
    #Log-Likelihood computation
    fbc=fb$c
    V=fb$beta * Em
    alpha=fb$alpha
    
    
    # tran_w =replicate(n_bin, matrix(0,n_states,n_states), simplify = FALSE)
    for(i  in 1:n) {
      tmp=fb$beta[i,1,]*alpha[i,1,]
      p=tmp/sum(tmp)
      ret[i,1]=niaidMI:::.samp(p)
      for(t in 1:(n_day-1)) {
        O1=alpha[i,t,]*Tran[[ bin[t] ]]
        tran_w0=niaidMI:::col_mult(O1, V[i,t+1,])
        p=tran_w0[ret[i,t],]
        p=p/sum(p)
        ret[i,t+1]=niaidMI:::.samp(p)
        #cat(globalcnt, i, t, p, "\n" )
      }
    }
    
    return(ret)
  }


bin <- rep(1, 27)

set.seed(2021)
#dataset_NM=sim_data(n=200)
dataset_NM=sim_data(n=200, drop_out_rate = 0.02, sporatic_rate = 0.2)

#Reformat dataset from dataset_NM to dataset_CH
dataset_CH <- NM2CH_data(dataset_NM)

start_BW <- get_start(dataset_CH[[1]][, -1], bin)
start_initP <- start_BW[[2]]
start_tP <- start_BW[[1]]

s=rep(1, nrow(dataset_NM))
Em = get_emission(dataset_NM, paste0("D",1:28))

# test1 <- .impute_one_debug(Pri=start_initP, s=s, Tran=start_tP, Em=Em, bin=bin)
# impute_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin = bin, start_initP, start_tP)

fit_CH <- BW_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin, start_initP, start_tP, message = FALSE, tol_llk = 1e-10)

initP <- fit_CH[[1]]
tP <- fit_CH[[2]]

NMimpute <- function(){
  NM <- NA
  for (i in 1:100) {
    NM <- rbind(NM, .impute_one_debug(Pri=list(initP), s=s, Tran=tP, Em=Em, bin=bin))
  }
  NM[-1, ]
}

CHimpute <- function(){
  CH <- NA
  for (i in 1:100) {
    CH <- rbind(CH, impute_CH(dataset_CH[[1]][, -1], dataset_CH[[2]], bin = bin, initP, tP))
  }
  CH[-1, ]
}

set.seed(12324)
test1 <- NMimpute()
set.seed(12324)
test2 <- CHimpute()


#table(test1 == test2)

test_that("Imputed values are all equal between BW and CH",
          {
            expect_equal(test1, test2)
          })


