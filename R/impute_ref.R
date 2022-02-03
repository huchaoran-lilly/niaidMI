#' @title Multiple Imputation for NIAID-OS using a reference.
#' @description Imputes NIAID OS data using a Markov model.
#' @seealso \code{\link{bootstrap_param_est}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). 
#' @param m Number of imputations.
#' @param ref Character vector with 1 column names. Specifies the reference group column. The reference group column must be logical, and there must be some (>=20) members of the reference group in each strata.
#' @param by Character vector with column names. Data will be broken up and imputed separately for every combination of values for specified columns in the data.
#' @param days Names of the columns that contain the score for each day. Columns should be in sequential order.
#' @param bin The assigned bin for pooling together information across transitions. Must be a numeric vector of length=(length(days)-1). By default all transitions are pooled together.
#' @param Em Emission probabilities. Generally the default should not be changed.
#' @param listFormatOut Return each imputed dataset in a list or combine into a single dataset.
#' @param tol Tolerance for relative reduction the log-likelihood to determine convergence of the Baum-Welch algorithm.
#' @param maxiter Maximum iterations before stopping the EM algorithm.
#' @param silent Allows silencing some messages.
#' @details
#' States for each patient/day in 'wide' may be the following: 
#' \itemize{
#'  \item{Not missing:}{An integer from 1 to 8.}
#'  \item{Missing:}{NA}
#'  \item{Partially Missing:}{ Range which may be code as a characters string such as '[1,7]' or '[1,2]'. Such a character string indicates that while the actual value is unknown, it is known that the value falls within the specified range. }
#' }
#' The reference based imputation uses a simple modification to the standard procedure.
#' First, within each strata, the model fit and bootstrapping procedure is performed only using the patients that are 
#' in the reference group. Second, the imputation for all of the patients in that strata (both reference and treated patients) 
#' is performed using the parameters as estimated/simulated based on the patients in the reference group.
#' 
#' @return 
#' If listFormatOut = TRUE, then a list will be returned with each element being an imputed data set.
#' If listFormatOut = FALSE, then a single data.frame will be returned where IMP_ID column is created.
#' 
#' 
#' @examples
#' test <- sim_data(100)
#' test$PBO=sample(c(TRUE, FALSE), size=nrow(test), replace = TRUE)
#' bs <- impute_ref(wide=test,ref="PBO",m=2, by="strata", silent=TRUE)
impute_ref <-
  function(wide, m, ref=NULL, by=NULL, days=paste0("D",1:28), bin=rep(1,length(days)-1), 
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
    
    if(is.null(ref))
      stop("ref must be specified.")
    if(!is.character(ref))
      stop("ref must be a character.")
    if(length(ref)!=1)
      stop("ref must be length 1.")
    if( !(ref %in% names(wide)) )
      stop("ref must be a column in the wide data.frame.")
    if( !is.logical(wide[[ref]])) 
      stop("ref must be a logical (i.e. TRUE/FALSE) column in the wide data.frame.")
    reftbl=table(factor(wide[[ref]],levels=c(TRUE,FALSE) ))
    if( reftbl["TRUE"]<20) 
      stop("All strata must have >=20 patients in the reference arm.")
    if( reftbl["FALSE"]==0) 
      stop("All patients are in the reference arm for at least one of the strata.")
    
    #originally I implemented to stratify baseline parameters. 
    #However, not planning to expose this functionality.
    #the following line sets up for no strata
    s=rep(1, nrow(wide))
    
    if(is.null(by)) { #no by grouping so just go
      ref_flg=wide[,ref]
      #Select only the reference data to fit the model
      wide_ref=wide[ref_flg,]
      fit=bootstrap_param_est(wide=wide_ref, b=m, bin=bin,
                              days=days, Em=Em[ref_flg,,],tol=tol,  
                              maxiter=maxiter, silent=silent)
      
      #impute all of the data using the reference model parameters
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
    
    if(fix_const!=0) {
      Pri=lapply(Pri,function(x) {
        x=x+fix_const
        x/sum(x)})
      Tran=lapply(Tran,function(x) {
        x=x+fix_const
        x/rowSums(x)})
    }
    
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

