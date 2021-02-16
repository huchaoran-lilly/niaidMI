### This file contains:
### fwd_bwd_CH: computes normalized forward/backward variables
### BW_CH     : estimates initial prob and transition prob
###             by Baum Welch algorithm
### impute_CH : imputes missing based on Markov chain (basic working block)
### imputeBS_CH:  imputes missing based on Markov chain and bootstrap
###               allowing strata imputation
### NM2CH     : transfer data from NM format to CH format
### get_start : get start value for BW algorithm according to NM's algorithm
### estTransMatrx : estimate transition matrix empirically
### estInitDist : estimate initial distribution empirically
### Check the examples at the bottom of the file for the
### explanation of arguments.

#######################################################################################################


## return a list with normalized forward/backward
## variables and c function
fwd_bwd_CH <- function(dat, q_func, bin, initP, tP){
  n_dat <- length(dat)
  alpha      <- matrix(NA, nrow = n_dat, ncol = 8) # not output/calculate
  alpha_star <- matrix(NA, nrow = n_dat, ncol = 8)
  beta       <- matrix(NA, nrow = n_dat, ncol = 8) # not output/calculate
  beta_star  <- matrix(NA, nrow = n_dat, ncol = 8)
  c_func     <- numeric(n_dat)
  
  alpha[1, ]         <- initP * q_func[dat[1], ]
  c_func[1]          <- sum(alpha[1, ])
  alpha_star[1, ]    <- alpha[1, ] / c_func[1]
  beta[n_dat, ]      <- 1
  beta_star[n_dat, ] <- 1
  
  for (i in 2:n_dat) {
    cart <- alpha_star[i - 1, ] %*% tP[[bin[i - 1]]] * q_func[dat[i], ]
    c_func[i] <- sum(cart)
    alpha_star[i, ] <- cart / c_func[i]
  }
  
  for (i in (n_dat - 1):1) {
    beta_star[i, ] <- tP[[bin[i]]] %*% (beta_star[i + 1, ] * q_func[dat[i + 1], ]) / c_func[i + 1]
  }
  
  return(list(alpha_star, beta_star, c_func))
}


#######################################################################################################


arraySum <- function(x) {
  result <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[2])
  for (i in 1:dim(x)[3]) {
    result <- result + x[, , i]
  }
  result
}


BW_CH <- function(dataset, q_func, bin, start_initP, start_tP,
                  tol_llk = 0, message = TRUE){
  n_data <- nrow(dataset)
  len_data <- ncol(dataset)
  n_bin <- length(unique(bin))
  initP  <- start_initP
  tP     <- start_tP
  ite <- 1
  llk_last <- 0
  
  while (TRUE) {
    ## initiate each iteration
    cart_updt_initP <- matrix(NA, nrow = n_data, ncol = 8)
    
    cart_updt_tP_nomi <- vector('list', n_bin)
    for (i in 1:n_bin) cart_updt_tP_nomi[[i]] <- array(NA, dim = c(8, 8, n_data))
    
    cart_updt_tP_deno <- vector('list', n_bin)
    for (i in 1:n_bin) cart_updt_tP_deno[[i]] <- matrix(NA, ncol = 8, nrow = n_data)
    
    llk <- 0

    for (i in 1:n_data) { ## calculate i-th patient
      ## get fwd bwd c_func
      cart <- fwd_bwd_CH(dataset[i, ], q_func, bin, initP, tP)
      
      ## get llk
      llk <- llk + sum(log(cart[[3]]))
      
      ## get gamma
      gamma <- cart[[1]] * cart[[2]]
      ## get xi
      ## test : each sum of each xi(t) is 1
      xi <- array(NA, dim = c(8, 8, len_data - 1))
      for (j in 1:(len_data - 1)) {
        cartA <- t(t(cart[[1]][j, ] * tP[[bin[j]]]) * cart[[2]][j + 1, ])
        cartA <- t(t(cartA) * q_func[dataset[i, j + 1], ]) / cart[[3]][j + 1]
        xi[, , j] <- cartA
      }
      
      # store result
      cart_updt_initP[i, ] <- gamma[1, ]
      
      for (s in 1:n_bin) {
        cart_updt_tP_nomi[[s]][, , i] <- arraySum(xi[, , which(bin == s)])
      }
      
      for (s in 1:n_bin) {
        cart_updt_tP_deno[[s]][i, ] <- apply(gamma[which(bin == s), ], 2, sum)
      }
    }
    
    ## update initP and transition prob
    initP <- apply(cart_updt_initP, 2, sum) / n_data
    for (s in 1:n_bin) {
      tP[[s]] <- arraySum(cart_updt_tP_nomi[[s]]) / apply(cart_updt_tP_deno[[s]], 2, sum)
    }
    
    ## fix state 8 related quant for NIAID
    #initP[8] <- 0
    for (s in 1:n_bin) tP[[s]][8, ] <- c(rep(0, 7), 1)
    
    ## print option
    if (message) print(paste0("Iteration ", ite, ": log-likelihood is ", llk, "."))
    
    ## check whether stop
    if (abs(llk_last - llk) <= tol_llk & ite >= 5) break
    
    ## update llk and iteration number
    llk_last <- llk
    ite <- ite + 1
  }
  
  return(list(initial_prob = initP,
              transition_prob = tP,
              iteration = ite,
              llk = llk))
}



#######################################################################################################


impute_CH <- function(dataset, q_func, bin, initP, tP){
  n_data <- nrow(dataset)
  len_data <- ncol(dataset)
  n_bin <- length(tP)
  
  imp_data <- matrix(NA, ncol = len_data, nrow = n_data)
  
  for (i in 1:n_data) { ## calculate i-th patient
    ## get fwd bwd c_func
    cart <- fwd_bwd_CH(dataset[i, ], q_func, bin, initP, tP)
    
    ## get gamma
    gamma <- cart[[1]] * cart[[2]]
    
    ## get xi
    ## test : each sum of each xi(t) is 1
    xi <- array(NA, dim = c(8, 8, len_data - 1))
    for (j in 1:(len_data - 1)) {
      cartA <- t(t(cart[[1]][j, ] * tP[[bin[j]]]) * cart[[2]][j + 1, ])
      cartA <- t(t(cartA) * q_func[dataset[i, j + 1], ]) / cart[[3]][j + 1]
      xi[, , j] <- cartA
    }
    
    ## impuate i-th patient
    imp_cart <- numeric(len_data)
    imp_cart[1] <- sample(1:8, size = 1, replace = TRUE, prob = gamma[1, ])
    for (j in 2:len_data) {
      prob_cart <- xi[imp_cart[j - 1], , j - 1] / sum(xi[imp_cart[j - 1], , j - 1])
      imp_cart[j] <- sample(1:8, size = 1, replace = TRUE, prob = prob_cart)
    }
    
    imp_data[i, ] <- imp_cart
  }
  
  imp_data
}



#####################################################################################################

## no strata
enc_imputeBS_CH <- function(data, q_func, bin, start_initP, start_tP, m,
                            tol_llk = 0, fix_const = 1E-6) { #m: num of bootstrap
  
  #fitMC <- BW_CH(data, q_func, bin, start_initP, start_tP, tol_llk = tol_llk, message = FALSE)
  
  result <- vector("list", m)
  cart <- result
  for (i in 1:m) {
    id <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
    cart[[i]] <- data[id, ]
  }
  
  for (i in 1:m) {
    fitMC <- BW_CH(cart[[i]], q_func, bin, start_initP, start_tP, tol_llk = tol_llk, message = FALSE)
    
    ## solve conflict of bootstrap and fwd/bwd calculation
    cart_initP <- fitMC[[1]] + fix_const
    cart_initP <- cart_initP / sum(cart_initP)
    
    cart_tP <- lapply(fitMC[[2]], function(x) {
      x <- x + fix_const
      x / rowSums(x)
    })
    ## this approach may be improved in future
    
    result[[i]] <- impute_CH(data, q_func, bin, cart_initP, cart_tP)
  }
  result
}


### Imputation function with bootstrap and Markov Chain model
## The procedure contains three steps: 1. fit Markov Chain with original
## dataset; 2. get bootstrap datasets from original dataset;
## 3. impute each bootstrap dataset separately with fitted Markov Chain.
## param:
##   data:     matrix contains original data. One column is strata if by != NULL.
##   q_func:   matrix with 0-1 specify q function.
##   bin:      vector inidcates bins cross time.
##   start_initP: start value of initial probability for BW algorithm to fit Markov chain.
##   start_tP: start value of transition probability for BW algorithm to fit Markov chain.
##   m:        int, number of imputations.
##   by:       int or char, indicates which column in data is strata.
##             If NULL, impute all together.
##   tol_llk:  num, indicates the tolerance of BW algorithm.
##   fix_const: num, a small enough number to solve conflict of bootstrap and real data.
## value:
##   a list with length m. Each component is one imputation.
## details:
##   If by != NULL, data must include one column indicating strata and column
##   name/num must be indicated in by. If by == NULL, then data must NOT contain
##   strata column.

## with strata
imputeBS_CH <- function(data, q_func, bin, start_initP, start_tP, m,
                        by = NULL, tol_llk = 0, fix_const = 1E-6) {
  if (is.null(by)) {
    result <- enc_imputeBS_CH(data, q_func, bin, start_initP, start_tP, m, tol_llk, fix_const)
    return(result)
  }
  
  result <- vector("list", m)
  
  by <- data[, by]
  unq_by <- unique(by)
  databy <- vector("list", length(unique(by)))
  for (i in seq_len(length(unique(by)))) {
    databy[[i]] <- data[by == unq_by[i], -1]
    cart <- enc_imputeBS_CH(databy[[i]], q_func, bin, start_initP, start_tP, m, tol_llk, fix_const)
    for (j in 1:m) {
      result[[j]] <- rbind(result[[j]], cbind(unq_by[i], cart[[j]]))
    }
  }
  
  result
}


#####################################################################################################


## 1-8 for 8 niaid score; 9 for NA; 10- for other
NM2CH_data <- function(data, days=paste0("D",1:28), strata = "strata") {
  ## 1-8 for 8 niaid score; 9 for NA; 10- for other
  dat <- data[, days]
  
  
  # get different states
  unq_state <- c()
  for (i in seq_len(ncol(dat))) {
    unq_state <- c(unique(as.character(dat[, i])), unq_state)
  }
  unq_state <- unique(unq_state)
  unq_state <- unq_state[!unq_state %in% as.character(1:8) & !is.na(unq_state)]
  unq_state_num <- seq_len(length(unq_state)) + 9
  
  
  # get q function
  q_func <- matrix(0, ncol = 8, nrow = length(unq_state) + 9)
  q_func[1:8, 1:8] <- diag(8)
  q_func[9, ] <- 1
  for (i in seq_len(length(unq_state))) {
    bg <- substr(unq_state[i], 2, 2)
    ed <- substr(unq_state[i], 4, 4)
    q_func[i + 9, seq(bg, ed)] <- 1
  }
  
  
  # get data
  tdata <- matrix(NA, ncol = length(days), nrow = nrow(data))
  for (i in seq_len(ncol(dat))) {
    for (j in seq_len(nrow(dat))) {
      
      cart_chr <- as.character(dat[j, i])
      if (cart_chr %in% as.character(1:8)) {
        tdata[j, i] <- as.numeric(cart_chr)
      } else {
        if (is.na(cart_chr)) {
          tdata[j, i] <- 9
        } else {
          tdata[j, i] <- unq_state_num[which(unq_state == cart_chr)]
        }
      }
    }
  }
  
  strt <- as.numeric(as.factor(data[, strata]))
  
  list(cbind(strt, tdata), q_func,
       cbind(unq_state, unq_state_num))
  
}



# ####################################################################################################

## get start value according to Nathan's code
get_start <- function(data, bin) {
  
  M = t(rbind(apply(data, 1, function(x) {
    x[x > 8] = NA
    x
  })))
  
  #originally I implemented to stratify baseline. However, not planning to expose this functionality.
  #the following line sets up for no strata
  s=rep(1, nrow(data))
  
  
  strt=niaidMI:::.get_start(M=M,s=s,bin=bin)
  
  list(strt$Tran, strt$Pri[[1]])
}

# ####################################################################################################


## the data matrix should not contains Subject ID
estTransMatrx <- function(data){ ## for niaid os only
  
  n_states <- length(unique(as.vector(data)))
  
  nomi <- matrix(0, ncol = n_states, nrow = n_states)
  deno <- matrix(0, ncol = n_states, nrow = n_states)
  
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data) - 1)) {
      if (is.na(data[i, j]) | is.na(data[i, j + 1])) break
      nomi[data[i, j], data[i, j + 1]] <- nomi[data[i, j], data[i, j + 1]] + 1
      deno[data[i, j], ] <- deno[data[i, j], ] + 1
    }
  }
  
  #print(nomi)
  #print(deno)
  nomi / deno
  
}

# estTransMatrx(imp_CH)
# estTransMatrx(imp_NM)
# 
# estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1])
# estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])
# 
# estTransMatrx(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1])
# estTransMatrx(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])




estInitDist <- function(data) { ## for niaid os only
  deno <- length(na.omit(data[, 1]))
  result <- numeric(8)
  for (i in 1:8) {
    result[i] <- sum(data[, 1] == i, na.rm = TRUE) / deno
  }
  result
}

# estInitDist(imp_CH)
# estInitDist(imp_NM)
# 
# estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 1), -1])
# estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 1), -1])
# 
# estInitDist(imp_strat_CH[which(imp_strat_CH[, 1] == 2), -1])
# estInitDist(imp_strat_NM[which(imp_strat_NM[, 1] == 2), -1])



# ####################################################################################################


################# EXAMPLE START HERE #################################################################
# ####################################################################################################
# 
# 
# NMdata <- niaidMI::sim_data(n=200)
# CHdata <- NM2CH_data(NMdata)
# 
# bin <- c(rep(1, 6), #(Please note data has 28 days)
#          rep(2, 7), rep(3, 7), rep(4, 7))
# 
# start_BW <- get_start(CHdata[[1]][, -1], bin)
# start_initP <- start_BW[[2]]
# start_tP <- start_BW[[1]]
# 
# BW_CH(CHdata[[1]][, -1], CHdata[[2]], bin, start_initP, start_tP)
# testimp  <- imputeBS_CH(CHdata[[1]], CHdata[[2]], bin, start_initP, start_tP, m = 3, by = 1)
# testimp2 <- imputeBS_CH(CHdata[[1]][, -1], CHdata[[2]], bin, start_initP, start_tP, m = 3)
#
#
# 
# ####################################################################################################
#
#
#
#
# ### example 1 (only 1 additional state 9 for missing)
# 
# ## q functions: a n_q_func by 8 matrix with (i, j)th
# ## element for q(h^obs=i, h=j)
# q_func <- rbind(diag(8), 1)
# 
# ## dat: observed NIAID for one patient
# dat <- sample(c(1:7, 9), size = 28, replace = TRUE)
# 
# ## bin: a vector with length(dat) - 1
# bin <- c(rep(1, 6), #(Please not data has 28 days)
#          rep(2, 7), rep(3, 7), rep(4, 7))
# 
# ## initial probability of chain
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
# 
# ## test fwd_bwd_CH
# test <- fwd_bwd_CH(dat, q_func, bin, initP, tP)
# apply(test[[1]]*test[[2]], 1, sum) ## should be 1
# 
# 
# ## dataset: observed NIAID for all patients
# ## each row is the record of one patient
# dataset <- matrix(NA, ncol = 28, nrow = 100)
# for (i in 1:100) {
#   dataset[i, ] <- sample(c(1:7, 9), size = 28, replace = TRUE)
# }
# dataset_fmt=dataset
# dataset_fmt[dataset_fmt==9]=NA
# colnames(dataset_fmt)=paste0("D",1:28)
# dataset_fmt=data.frame(dataset_fmt)
# 
# 
# ## a quick example
# result1 <- BW_CH(dataset, q_func, bin, start_initP, start_tP)
# result2 = bootstrap_param_est(wide=dataset_fmt, b=2, bin=bin, tol = 1e-20)[[1]]
# result1$initial_prob
# result2$Pri
# 
# result1[[2]]
# result2$Pri
# 
# 
# 
# ## quick example
# impData <- impute_CH(dataset, q_func, bin, result1[[1]], result1[[2]])
#
# ##########################################################################
# 
# 
# ### example 2 (2 additional states, 9 for missing, 10 for missing not die)
# 
# ## q functions: a n_q_func by 8 matrix with (i, j)th
# ## element for q(h^obs=i, h=j)
# q_func <- rbind(diag(8), 1, c(rep(1, 7), 0))
# 
# ## dat: observed NIAID for one patient
# dat <- sample(c(1:7, 9, 10), size = 28, replace = TRUE)
# 
# ## bin: a vector with length(dat) - 1
# bin <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7))
# 
# ## initial probability of chain
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
# 
# ## test fwd_bwd_CH
# test <- fwd_bwd_CH(dat, q_func, bin, initP, tP)
# apply(test[[1]]*test[[2]], 1, sum) ## should be 1
# 
# 
# ## dataset: observed NIAID for all patients
# ## each row is the record of one patient
# dataset <- matrix(NA, ncol = 28, nrow = 100)
# for (i in 1:100) {
#   dataset[i, ] <- sample(c(1:7, 9, 10), size = 28, replace = TRUE)
# }
# 
# 
# ## a quick example
# result <- BW_CH(dataset, q_func, bin, start_initP, start_tP)
# 
# ## quick example
# impData <- impute_CH(dataset, q_func, bin, result[[1]], result[[2]])
# 
# ## should be all true
# (dataset %in% c(9, 10)) != (impData == dataset)
# 
# ##########################################################################


