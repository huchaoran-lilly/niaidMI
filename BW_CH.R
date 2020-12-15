## q functions: a n_q_func by 8 matrix with (i, j)th
## element for q(h^obs=i, h=j)
q_func <- rbind(diag(8), 1)

## dat: observed NIAID for one patient
dat <- sample(c(1:7, 9), size = 29, replace = TRUE)

## bin: a vector with length(dat) - 1
bin <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7))

## initial probability of chain
initP <- c(rep(1, 7) / 7, 0)

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


test <- fwd_bwd_CH(dat, q_func, bin, initP, tP)
apply(test[[1]]*test[[2]], 1, sum)



















