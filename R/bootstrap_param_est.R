#' @title Estimation of Markov model.
#' @description Fits a Markov model then bootstraps the data and refits the model.
#' @seealso \code{\link{impute}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). See details.
#' @param b Number of bootstrap samples to take.
#' @param days Names of the columns that contain the score for each day.
#' @param bin The assigned bin for pooling together information across transitions. Must be a numeric vector of length=(length(days)-1). By default all transitions are pooled together.
#' @param Em Emission probabilities. Default should be used unless user is advanced.
#' @param tol Tolerance for relative reduction the log-likelihood to determine convergence of the Baum-Welch algorythm.
#' @param maxiter Maximum iterations before stopping the EM algorithm.
#' @param silent Allows silencing some messages.
#' @details
#' States for each patient/day in 'wide' may be the following: 
#' \itemize{
#'  \item{Not missing:}{An integer from 1 to 8.}
#'  \item{Missing:}{NA}
#'  \item{Partially Missing:}{ range which may be code as a characters string such as '[1,7]' or '[1,2]'. Such a character string indicates that while the actual value is unknown, it is known that the value falls within the specified range. }
#' }
#' Generally the user will not need to call this function directly because it is called by the 'impute' function.
#' 
#' @return 
#' A list object with the following components:
#' \describe{
#'   \item{fit}{Contains results of the primary model fit}
#'   \item{boot}{Contains relults from the bootstrap model fit.}
#'   \item{bin}{The input.}
#'   \item{s}{Ignor. May be used in future.}
#'   \item{days}{From input.}
#'   \item{Em}{From input.}
#' }
#' 
#' 
#' @examples
#' test <- sim_data(100)
#' bs <- bootstrap_param_est(wide=test,b=2)

bootstrap_param_est <-
function(wide, b, days=paste0("D",1:28), bin=rep(1,length(days)-1), 
         Em=get_emission(wide, days), tol=1E-6, maxiter=200, silent=FALSE) {

  
  if(!is.data.frame(wide))
    stop("wide must be a data.frame.")

  if(!is.numeric(bin))
    stop("bin must be numeric.")
  if(!is.vector(bin))
    stop("bin must be a vector")
  if(length(bin)!=(length(days)-1))
    stop("length(bin) must be the same as length(days)-1.")
  
  if(!is.numeric(b))
    stop("b must be numeric.")
  if(!is.vector(b))
    stop("b must be a vector")
  if(length(b)!=1)
    stop("length(b) must be 1.")
  if(b<0)
    stop("b must be >=0.")
  
  M = wide[,days]
  M = as.matrix(do.call(cbind, lapply(M, function(x) {
    x=suppressWarnings(as.numeric(as.character(x)))
    # x[x<1 | x>8]=NA
    x
  })))
  
  #originally I implemented to stratify baseline. However, not planning to expose this functionality.
  #the following line sets up for no strata
  s=rep(1, nrow(wide))
  
  
  strt=.get_start(M=M,s=s,bin=bin)
  if(!silent) cat("Fiting Model\n")
  fit=.baum_welsh(Pri=strt$Pri, s=s, Tran=strt$Tran, Em=Em, bin=bin, tol=tol,  maxiter=maxiter)
  
  ret=list(fit=fit, s=s, bin=bin, days=days, Em=Em)
  
  if(!silent) cat("Fiting Bootstrap Samples:\n")
  ret$boot=lapply(1:b, function(i) {
    if(!silent) cat("    b=",i)
    boot=.bootstrap_dta(Em=Em, M=M, s=s)
    # strt=.get_start(M=boot$M,s=boot$s, bin=bin)
    
    fit1=.baum_welsh(Pri=strt$Pri, s=boot$s, Tran=fit$Tran, Em=boot$Em, bin=bin, tol=tol,  maxiter=maxiter)
    # fit1=.baum_welsh(Pri=strt$Pri, s=boot$s, Tran=fit$Tran, Em=boot$Em, bin=bin, tol=tol,  maxiter=maxiter)
    if(!silent) cat(", iterations=", fit1$nit,"\n")
    return(fit1[c("Tran","Pri")])
  })
  
  return(ret)
}

.bootstrap_dta <-
  function(Em, M, s) {
    pos=sample(1:nrow(M), replace = TRUE)
    #cat(pos, "\n")
    return(list(Em=Em[pos,,], M=M[pos,], s=s[pos]))
  }
.samp <-
  function(p) {
    sample.int(8,1,replace=TRUE,prob=p)
  }


.get_start <-
function(M,s,bin) {

  Tran=lapply(split(1:length(bin), bin), function(i) {
    .get_tran(M, min(i):(max(i)+1))
  })

  Pri=lapply(split(M[,1],s), function(x) {
    ret=table(factor(x,levels=1:8))
    as.vector(ret/sum(ret))
  }
  )
  
  return(list(Tran=Tran, Pri=Pri))
}


.get_tran <-
function(M, pos, const=1E-6) {

  M=M[,pos]
  tbl=lapply(1:nrow(M), function(i) {
    v=M[i,]
    x2=factor(v[-1], levels=1:8)
    x1=factor(v[-length(v)], levels=1:8)
    tbl=table(x1, x2)
    if(!all(tbl[8,1:7]==c(0,0,0,0,0,0,0)))
      stop("Data error. Patient transitioned from state 8 (dead) to another state. No resurections allowed in the data.")
    return(tbl)
  })
  tbltot=tbl[[1]]
  for(i in 2:length(tbl))
    tbltot=tbltot+tbl[[i]]

  tranNum=as.matrix(tbltot) + const 
  tranNum=tranNum/rowSums(tranNum)
  tranNum[8,]=c(0,0,0,0,0,0,0,1) #don't apply constant to 8
  
  return(tranNum)

}



