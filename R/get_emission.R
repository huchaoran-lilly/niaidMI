#' @title Emission probability evaluation.
#' @description Creates emission probabilities from data. Generally for internal niaidMI package use only.
#' @seealso \code{\link{impute}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). 
#' @param days Names of the columns that contain the score for each day.
#' @details
#' States for each patient/day in 'wide' may be the following: 
#' \itemize{
#'  \item{Not missing:}{An integer from 1 to 8.}
#'  \item{Missing:}{NA}
#'  \item{Partially Missing:}{ Range which may be code as a characters string such as '[1,7]' or '[1,2]'. Such a character string indicates that while the actual value is unknown, it is known that the value falls within the specified range. }
#' }
#' Generally the user will not need to call this function directly because it is called by the 'impute' function.
#
#' @return 
#' Creates a 3 dimensional array that is "number of patients" x "number
#' of days" x "8 NIAID stats." This array contains only 1 or 0 for each entry indicating 
#' if the state for a given day and individual is consistent with the data.
#' 
#' @examples
#' test <- sim_data(200)
#' Em <- get_emission(wide=test,days=paste0("D",1:28))

get_emission=function(wide,days) {
  if(!is.data.frame(wide))
    stop("wide must be a data.frame.")
  
  tmp=.check_bad_transition(wide,days)
  
  split_list=list()
  for(nm in days) {
    col_str=wide[[nm]]
    splt=sub("[", "", col_str, fixed=TRUE)
    splt=sub("]", "", splt, fixed=TRUE)
    splt=strsplit( splt, ",", fixed=TRUE)
    split_list[[nm]]=mapply(function(x,v) {
      
      ret=tryCatch(as.numeric(x), 
                   warning=function(e) stop(paste("Could not process niad os value:", v,".")))
      if(any(is.na(ret))) {
        if(length(ret)>1)
          stop(paste("Could not process niad os value of:", v,"."))
        return(c(1,8)) #could be anything
      } else if(any(ret<1 | ret>8 | round(ret)!=ret)){
        stop(paste("Niad os value must be integer between 1 and 8. Found value of:", v,"."))
      } else if(length(ret)>2) {
        stop(paste("Could not process niad os value of:", v,"."))
      } else if(length(ret)==2) {
        if(ret[1]==ret[2])
          stop(paste("Could not process niad os value of:", v,"."))
        return(c(min(ret),max(ret)))
      } else
        return(ret)
    }, splt,col_str, SIMPLIFY = FALSE)
  }
  
  valuelist=unique(do.call(c, split_list))
  
  tochar=function(x) sapply(x, function(x) paste(x,collapse="_"))
  
  emmssionNum=lapply(valuelist, function(x) {
    ret=rep(0,8)
    if(length(x)==1) {
      ret[x]=1
    } else {
      ret[x[1]:x[2]]=1
    }
    return(ret)
  })
  names(emmssionNum)=tochar(valuelist)
  Em=array(0,c(nrow(wide),length(days),8))
  Emlst=lapply(split_list, function(x) do.call(rbind,emmssionNum[tochar(x)]))
  for(i in 1:length(days))
    Em[,i,]=Emlst[[i]]
  return(Em)
}

.check_bad_transition=function(wide,days) {
  tmp=wide[,days]
  for(i in 1:nrow(wide)) {
    v = unlist(tmp[i,])
    pos8 = which(v==8)
    if(length(pos8)>0) {
      first=pos8[1]
      flgmaybe8=grepl("8", v, fixed=TRUE) | is.na(v)
      if(any(!flgmaybe8[first:length(days)]))
        stop(paste("Cannot transition out of 8:\n", 
                   paste(wide[i,], collapse=",")))
    }
  }
  return(TRUE)
}
