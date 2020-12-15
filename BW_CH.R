## q functions: a n_q_func by 8 matrix with (i, j)th
## element for q(h^obs=i, h=j)
q_func <- rbind(diag(8), 1)

## dat: observed NIAID for one patient
dat <- sample(1:9, size = 29, replace = TRUE)

## bin: a vector with length(dat) - 1
bin <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7))

## initial probability of chain
init <- c(rep(1, 7) / 7, 0)

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


## return a matrix with 6 columns for
## time, alpha, alpha*, c, beta, beta*
fwd_bwd_CH <- function(dat){
  result <- matrix(NA, nrow = length(dat), ncol = 6)
}

























