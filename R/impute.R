#' @title Function impute
#' @description Imputes NIAID OS data using a hidden markov model.
#' @seealso \code{\link{bootstrap_param_est}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). 
#' Data must be in extended representation. Here as state of (9) is missing and possibly dead, and 
#' a state of (10) is missing but not dead.
#' @param bin The assigned bin for each transition. Must be a numeric vector of length=(1-length(days)).
#' @param m Number of imputations.
#' @param by Character vector with column names. Data will be broken up and imputed separately for every combination of values for those columns in the data.
#' @param baselineStrat Character vector with column names. Baseline proportion estimated separately for each strata.
#' @param Em Emission probabilities. Generally the default should not be changed.
#' @param days Names of the columns that contain the score for each day.
#' @param wideFormatOut Return in wide format (TRUE) or long format (FALSE)
#' @param listFormatOut Return each imputed dataset in a list or combine into a single dataset.
#' @param tol tolerance for relative reduction the log-likelihood to determine convergence of the Baum-Welch algorythm.
#' @param silent Allows silencing some messages.
#' @param s For internal use only.
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
impute <-
function(wide, bin, m, by=NULL, baselineStrat=NULL, days=paste0("D",1:28),
                  Em=get_emission(wide, days), wideFormatOut=TRUE,
                  listFormatOut=TRUE, tol=1E-6, silent=FALSE, s=NULL) {
  
  #Setup up baseline starta assigning each patient to a strata
  if(is.null(baselineStrat)) {
    s=rep(1, nrow(wide))
  } else {
    if(!is.null(s))
      stop("Cannot specify both 'baselineStrat' and 's.'")
    if(!is.character(baselineStrat))
      stop("'baselineStrat' must be of type character.")
    if(!all(baselineStrat %in% colnames(wide)))
      stop("All variables listed in by must be in the 'wide' data.frame.")
    if(is.null(by))
      if(any(baselineStrat %in% by))
        stop("Baseline stratification factors cannot be included as a 'by' grouping.")
    tmp=by_str_for_impute=apply(wide[,baselineStrat], 1, function(X) paste(X,collapse="|"))
    s=as.numeric(as.factor(tmp))
  }

  
  if(is.null(by)) { #no by grouping so just go
    fit=bootstrap_param_est(wide=wide,  s=s, bin=bin, b=m,
                            days=days, Em=Em,tol=tol, silent=silent)
    
    imp=lapply(fit$boot, function(bfit) {
      .impute_one(Pri=bfit$Pri, s=s, Tran=bfit$Tran, 
                 Em=Em, 
                 bin=bin)
    })
    extravars=wide[,which(!(colnames(wide) %in% days))]
    if(wideFormatOut) {
      ret=mapply(function(impi,i) {
        colnames(impi)=days
        data.frame(IMP_ID=i,extravars,impi)
      }, imp, 1:length(imp), SIMPLIFY = FALSE)
    } else {
      orig=melt(t(as.matrix(wide[,days])))
      ret=mapply(function(impi,i) {  
        colnames(impi)=days
        rownames(impi)=1:nrow(impi)
        impi=melt(t(impi))
        names(impi)=c("DAY", "PATIENT", "value")
        impi$orig_value=orig[,3]
        data.frame(IMP_ID=i,extravars[imp$PATIENT,],imp)
      }, imp, 1:length(imp), SIMPLIFY = FALSE)
    }
  } else {

    if(!is.character(by))
      stop("'by' must be of type character.")
    if(!all(by %in% colnames(wide)))
      stop("All variables listed in by must be in the 'wide' data.frame.")
    
    by_str_for_impute=apply(wide[,by], 1, function(X) paste(X,collapse="|"))

    poslst=split(1:nrow(wide),by_str_for_impute)
    imp=lapply(poslst, function(pos) {
      # browser()
      wide_by=wide[pos,]
      s_by=s[pos]
      Em_by=Em[pos,,]
      impute(wide=wide_by,  s=s_by, bin=bin, m=m, 
             days=days, Em=Em_by, 
             wideFormatOut=wideFormatOut, listFormatOut=TRUE, 
             by=NULL, silent=silent)
    })

    ret=list()
    for(i in 1:m)
      ret[[i]]=do.call(rbind, lapply(imp,function(x) x[[i]] ))
    
    ret=lapply(ret, function(x) {
      rownames(x)=NULL
      return(x)
    })
  }
  if(!listFormatOut) {
    ret=do.call(rbind, ret)
  }
  return(ret)
}




.impute_one <-
  function(Pri, s, Tran, Em, bin, n_bin=max(bin),
           n=dim(Em)[1], n_day=dim(Em)[2], n_states=8, n_strata=max(s), fix_const=1E-6) {
    
    
    Pri=lapply(Pri,function(x) {
      x=x+fix_const
      x/sum(x)})
    Tran=lapply(Tran,function(x) {
      x=x+fix_const
      x/rowSums(x)})
    
    fb <- .forward_backward(Pri=Pri, s=s, Tran=Tran, Em=Em, bin=bin,
                           n_day=n_day, n_states=n_states, n_bin=n_bin, n_strata=n_strata)
    
    ret=matrix(0, n, n_day)
    for(i in 1:n) {
      
      # cat(i, "\n")
      #day 1
      tmp=fb$alpha[i,1,]*fb$beta[i,1,]
      p=tmp/sum(tmp)
      ret[i,1]=.samp(p)
      
      #other days
      for(t in 1:(n_day-1)) {
        # cat(" ", t, "\n")
        O1=fb$alpha[i,t,]*Tran[[ bin[t] ]]
        v1=fb$beta[i,t+1,] * Em[i,t+1,]
        tran_w0=col_mult(O1, v1)
        p=tran_w0[ret[i,t],]
        p=p/sum(p)
        ret[i,t+1]=.samp(p)
      }
    }
    return(ret)
  }

