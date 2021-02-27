#' @title Simulate data.
#' @description Simulate data for the purpose of testing the package.
#' @export
#' @param n Number of samples. 
#' @param fit Contains transition and prior parameters.
#' @param strata Named vector with probabilities to simulate strata.
#' @param days Names of the output columns that contain the score for each day.
#' @param bin The assigned bin for pooling together information across transitions. Must be a numeric vector of length=(length(days)-1). By defualt all transitions are pooled together.
#' @param drop_out_rate Parameter for simulating dropout. Rate is per day.
#' @param sporatic_rate Parameter for simulating missing.
#' @param miss_not_dead_rate Parameter to control missing but not dead rate.
#' @details
#' This simulation function is used to generated data for examples and testing of the package.
#' 
#' @return 
#' Returns wide format data. Possible niad os values may be:
#' \itemize{
#'  \item{Not missing:}{An integer from 1 to 8.}
#'  \item{Missing:}{NA}
#'  \item{Partially Missing:}{ range which may be code as a characters string such as '[1,7]' or '[1,2]'. Such a character string indicates that while the actual value is unknown, it is known that the value falls within the specified range. }
#' }
#' @examples
#' test <- sim_data(200)
#' 

sim_data <-
  function(n, 
           fit=list(Pri=list(c(0,0,0,.5,.25,.25,0,0)),
                    Tran=list(matrix(c(0.74, 0.16, 0.04, 0.01, 0, 0, 0, 0, 0.19, 0.63, 0.15, 
                                       0.04, 0.01, 0, 0, 0, 0.05, 0.16, 0.61, 0.15, 0.04, 0.01, 0, 0, 
                                       0.02, 0.04, 0.15, 0.6, 0.15, 0.04, 0.01, 0, 0, 0.01, 0.04, 0.15, 
                                       0.6, 0.15, 0.04, 0, 0, 0, 0.01, 0.04, 0.15, 0.61, 0.16, 0, 0, 
                                       0, 0, 0.01, 0.04, 0.15, 0.63, 0, 0, 0, 0, 0, 0.01, 0.04, 0.16, 
                                       1), 8,8))),
           strata=c("s1"=.2, "s2"=.8),
           days=paste0("D",1:28), 
           bin=rep(1,length(days)-1),
           drop_out_rate=.01,
           sporatic_rate=.05,
           miss_not_dead_rate=.2
           
  ) {
    if(!is.numeric(bin))
      stop("bin must be numeric.")
    if(!is.vector(bin))
      stop("bin must be a vector")
    if(length(bin)!=(length(days)-1))
      stop("length(bin) must be the same as length(days)-1.")
    
    ret=matrix(0,n,length(days))
    colnames(ret)=days
    for(i in 1:n) {
      ret[i,1]=.samp(fit$Pri[[1 ]])
      for(t in 2:length(days)) {
        mytran=fit$Tran[[ bin[t-1] ]]
        ret[i,t]=.samp(mytran[ret[i,t-1],])
      }
    }
    mode(ret)="character"
    
    for(i in 1:n) {
      for(t in 1:length(days)) {
        if(!is.na(ret[i,t]) & 
           ret[i,t]!="[1,7]" &
           ret[i,t]!="8") {

            if(runif(1)<drop_out_rate){
              if(runif(1)>miss_not_dead_rate){
                for(t1 in t:length(days))
                  ret[i,t1]=NA
              } else {
                for(t1 in t:length(days))
                  ret[i,t1]="[1,7]"
              }
            }
          
          
          if(!is.na(ret[i,t]) & ret[i,t]!="[1,7]") {
            if(runif(1)<sporatic_rate){
              if(runif(1)>miss_not_dead_rate)
                ret[i,t]=NA
              else
                ret[i,t]="[1,7]"
            }
          }
        }
      }
      
    }
    ret=data.frame(ID=NA, 
               strata=sample(names(strata), n, prob=strata, replace = TRUE),
               ret)
    ret=ret[order(ret$strata, rowSums(is.na(ret[,days]) | ret[,days]=="[1,7]")),]
    rownames(ret)=NULL
    ret$ID=1:n
    return(ret)
  }



