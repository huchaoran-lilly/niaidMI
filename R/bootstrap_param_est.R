#' @title Function bootstrap_param_est
#' @description Fits a hidden markov model then bootstraps the data and refits the model.
#' @seealso \code{\link{tweights}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). 
#' Data must be in extended representation. Here as state of (9) is missing and possibly dead, and 
#' a state of (10) is missing but not dead
#' @param s Numeric vecor representing the assigned strata for each patient as a number.
#' @param bin The assigned bin for each transition. Must be a numeric vector of length=(1-length(days)).
#' @param days Names of the columns that contain the score for each day.
#' @param Em Emission probabilities.
#' @param tol tolerance for relative reduction the log-likelihood to determine convergence of the Baum-Welch algorythm.
#' @param silent Allows silencing some messages.
#' @details
#' 
#' @return 
#' A list object with the following components
#' \describe{
#'   \item{fit}{Contains results of the primary model fit}
#'   \item{boot}{Contains relults from the bootstrap model fit.}
#'   \item{bin}{The input.}
#'   \item{days}{Formatted dataset.}
#'   \item{Em}{Attempted target.}
#' }
#' 
#' 
#' @examples


bootstrap_param_est <-
function(wide,  s, bin, b, days=paste0("D",1:28), 
         Em=get_emission(wide, days), tol=1E-6, silent=FALSE) {

  
  M = wide[,days]
  M = as.matrix(do.call(cbind, lapply(M, function(x) {
    x[x>8]=NA
    x
  })))
  
  
  strt=.get_start(M=M,s=s,bin=bin)
  if(!silent) cat("Fiting Model\n")
  fit=.baum_welsh(Pri=strt$Pri, s=s, Tran=strt$Tran, Em=Em, bin=bin, tol=tol)
  if(!silent) cat("Fiting Bootstrap Samples:\n")
  boot=lapply(1:b, function(i) {
    if(!silent) cat("    b=",i,"\n")
    boot=.bootstrap_dta(Em=Em, M=M, s=s)
    strt=.get_start(M=boot$M,s=boot$s, bin=bin)
    fit=.baum_welsh(Pri=strt$Pri, s=boot$s, Tran=strt$Tran, Em=boot$Em, bin=bin, tol=tol)
    return(fit[c("Tran","Pri")])
  })
  
  return(list(fit=fit, boot=boot, s=s, bin=bin, days=days, Em=Em))
}

.bootstrap_dta <-
  function(Em, M, s) {
    pos=sample(1:nrow(M), replace = TRUE)
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
function(M, pos, const=1E-20) {

  M=M[,pos]
  tbl=lapply(1:nrow(M), function(i) {
    v=M[i,]
    x2=factor(v[-1], levels=1:8)
    x1=factor(v[-length(v)], levels=1:8)
    tbl=table(x1, x2)
  })
  tbltot=tbl[[1]]
  for(i in 2:length(tbl))
    tbltot=tbltot+tbl[[i]]

  tranNum=as.matrix(tbltot) + const 
  tranNum=tranNum/rowSums(tranNum)

  return(tranNum)

}




