#' @title Function impute
#' @description Imputes NIAID OS data using a hidden markov model.
#' @seealso \code{\link{bootstrap_param_est}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). 
#' @param m Number of imputations.
#' @param by Character vector with column names. Data will be broken up and imputed separately for every combination of values for those columns in the data.
#' @param days Names of the columns that contain the score for each day. Columns should be in sequential order.
#' @param bin The assigned bin for pooling together information across transitions. Must be a numeric vector of length=(length(days)-1). By defualt all transitions are pooled together.
#' @param Em Emission probabilities. Generally the default should not be changed.
#' @param listFormatOut Return each imputed dataset in a list or combine into a single dataset.
#' @param tol tolerance for relative reduction the log-likelihood to determine convergence of the Baum-Welch algorythm.
#' @param maxiter maximum iterations before stopping the EM algorithm.
#' @param silent Allows silencing some messages.
#' @details
#' States for each patient/day in 'wide' may be the following: 
#' \itemize{
#'  \item{Not missing:}{An integer from 1 to 8.}
#'  \item{Missing:}{NA}
#'  \item{Partially Missing:}{ range which may be code as a characters string such as '[1,7]' or '[1,2]'. Such a character string indicates that while the actual value is unknown, it is known that the value falls within the specified range. }
#' }
#' @return 
#' If listFormatOut = TRUE, then a list will be returned with each element being an imputed data set.
#' If listFormatOut = FALSE, then a single data.frame will be returned where IMP_ID column is created.
#' 
#' 
#' @examples
#' test <- sim_data(200)
#' bs <- impute(wide=test,m=5, by="strata", silent=TRUE)
impute <-
  function(wide, m, by=NULL, days=paste0("D",1:28), bin=rep(1,length(days)-1), 
           Em=get_emission(wide, days),
           listFormatOut=FALSE, tol=1E-6, maxiter=200, silent=FALSE) {

    
    if(!is.data.frame(wide))
      stop("wide must be a data.frame.")
    
    if(!is.numeric(bin))
      stop("bin must be numeric.")
    if(!is.vector(bin))
      stop("bin must be a vector")
    if(length(bin)!=(length(days)-1))
      stop("length(bin) must be the same as length(days)-1.")
    
    if(!is.numeric(m))
      stop("m must be numeric.")
    if(!is.vector(m))
      stop("m must be a vector")
    if(length(m)!=1)
      stop("length(m) must be 1.")
    if(m<2)
      stop("m must be >=2.")
    
    #originally I implemented to stratify baseline. However, not planning to expose this functionality.
    #the following line sets up for no strata
    s=rep(1, nrow(wide))
    
    
    if(is.null(by)) { #no by grouping so just go
      fit=bootstrap_param_est(wide=wide, b=m, bin=bin,
                              days=days, Em=Em,tol=tol,  
                              maxiter=maxiter, silent=silent)
      
      imp=lapply(fit$boot, function(bfit) {
        .impute_one(Pri=bfit$Pri, s=s, Tran=bfit$Tran, 
                    Em=Em, 
                    bin=bin)
      })
      extravars=wide[,which(!(colnames(wide) %in% days))]
      
      ret=mapply(function(impi,i) {
        colnames(impi)=days
        data.frame(IMP_ID=i,extravars,impi)
      }, imp, 1:length(imp), SIMPLIFY = FALSE)
      
    } else {
      
      if(!is.character(by))
        stop("'by' must be of type character.")
      if(!all(by %in% colnames(wide)))
        stop("All variables listed in by must be in the 'wide' data.frame.")
      
      by_str_for_impute=apply(wide[,by, drop=FALSE], 1, function(X) paste(X,collapse="|"))
      
      poslst=split(1:nrow(wide),by_str_for_impute)
      imp=lapply(poslst, function(pos) {
        # browser()
        wide_by=wide[pos,]
        s_by=s[pos]
        Em_by=Em[pos,,]
        impute(wide=wide_by, bin=bin, m=m, 
               days=days, Em=Em_by, 
               listFormatOut=TRUE, 
               by=NULL, tol=tol,  
               maxiter=maxiter,
               silent=silent)
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
    
    
    #Log-Likelihood computation
    fbc=fb$c
    V=fb$beta * Em
    alpha=fb$alpha


    # tran_w =replicate(n_bin, matrix(0,n_states,n_states), simplify = FALSE)
    for(i  in 1:n) {
      tmp=fb$beta[i,1,]*alpha[i,1,]
      p=tmp/sum(tmp)
      ret[i,1]=.samp(p)
      for(t in 1:(n_day-1)) {
        O1=alpha[i,t,]*Tran[[ bin[t] ]]
        tran_w0=col_mult(O1, V[i,t+1,])
        p=tran_w0[ret[i,t],]
        p=p/sum(p)
        ret[i,t+1]=.samp(p)
      }
    }
    
    return(ret)
  }

