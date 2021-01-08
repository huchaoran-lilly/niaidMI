### This file contains:
### fwd_bwd_CH: computes normalized forward/backward variables
### BW_CH     : estimates initial prob and transition prob
###             by Baum Welch algorithm
### impute_CH : imputes missing based on estimated process
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
    initP[8] <- 0
    for (s in 1:n_bin) tP[[s]][8, ] <- c(rep(0, 7), 1)
    
    ## print option
    if (message) print(paste0("Iteration ", ite, ": log-likelihood is ", llk, "."))
    
    ## check whether stop
    if (abs(llk_last - llk) <= tol_llk & ite >= 2) break
    
    ## update llk and iteration number
    llk_last <- llk
    ite <- ite + 1
  }
  
  return(list(initial_prob = initP,
              transition_prob = tP,
              iteration = ite))
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


#######################################################################################################




### example 1 (only 1 additional state 9 for missing)

## q functions: a n_q_func by 8 matrix with (i, j)th
## element for q(h^obs=i, h=j)
q_func <- rbind(diag(8), 1)

## dat: observed NIAID for one patient
dat <- sample(c(1:7, 9), size = 29, replace = TRUE)

## bin: a vector with length(dat) - 1
bin <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7))

## initial probability of chain
initP <- c(rep(1, 7) / 7, 0)
start_initP <- initP

## transition probability: a list with length equals
## to the number of bins. Each element is a 8*8 matrix.
tP <- vector('list', length(unique(bin)))
tP[[1]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
tP[[2]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
tP[[3]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
tP[[4]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
start_tP <- tP

## test fwd_bwd_CH
test <- fwd_bwd_CH(dat, q_func, bin, initP, tP)
apply(test[[1]]*test[[2]], 1, sum) ## should be 1


## dataset: observed NIAID for all patients
## each row is the record of one patient
dataset <- matrix(NA, ncol = 29, nrow = 100)
for (i in 1:100) {
  dataset[i, ] <- sample(c(1:7, 9), size = 29, replace = TRUE)
}


## a quick example
result <- BW_CH(dataset, q_func, bin, start_initP, start_tP)

## quick example
impData <- impute_CH(dataset, q_func, bin, result[[1]], result[[2]])


##########################################################################


### example 2 (2 additional states, 9 for missing, 10 for missing not die)

## q functions: a n_q_func by 8 matrix with (i, j)th
## element for q(h^obs=i, h=j)
q_func <- rbind(diag(8), 1, c(rep(1, 7), 0))

## dat: observed NIAID for one patient
dat <- sample(c(1:7, 9, 10), size = 29, replace = TRUE)

## bin: a vector with length(dat) - 1
bin <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7))

## initial probability of chain
initP <- c(rep(1, 7) / 7, 0)
start_initP <- initP

## transition probability: a list with length equals
## to the number of bins. Each element is a 8*8 matrix.
tP <- vector('list', length(unique(bin)))
tP[[1]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
tP[[2]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
tP[[3]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
tP[[4]] <- rbind(matrix(1/8, nrow = 7, ncol = 8),
                 c(rep(0, 7), 1))
start_tP <- tP

## test fwd_bwd_CH
test <- fwd_bwd_CH(dat, q_func, bin, initP, tP)
apply(test[[1]]*test[[2]], 1, sum) ## should be 1


## dataset: observed NIAID for all patients
## each row is the record of one patient
dataset <- matrix(NA, ncol = 29, nrow = 100)
for (i in 1:100) {
  dataset[i, ] <- sample(c(1:7, 9, 10), size = 29, replace = TRUE)
}


## a quick example
result <- BW_CH(dataset, q_func, bin, start_initP, start_tP)

## quick example
impData <- impute_CH(dataset, q_func, bin, result[[1]], result[[2]])

## should be all true
(dataset %in% c(9, 10)) != (impData == dataset)

##########################################################################


