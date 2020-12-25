#' @title Function get_emission
#' @description Creates Emission probabilities from data
#' @seealso \code{\link{impute}}
#' @export
#' @param wide Data in wide format (i.e., each day is a column). 
#' Data must be in extended representation. Here as state of (9) is missing and possibly dead, and 
#' a state of (10) is missing but not dead
#' @param days Names of the columns that contain the score for each day.
#' @details
#' Documentation of how this matrix is created is described in the vignette.
#' @return 
#' Creates a 3 dimensional array that is "number of patients" x "number
#' of days" x "8 NIAID stats." 
#' 
#' @examples
#' 
#' 


get_emission <-
function(wide, days) {
  #Setup Emission Values
  #  + 9 is missing and possibly dead 
  #  + 10 is missing but not dead
  emmssionNum=cbind(diag(8),rep(1,8),c(rep(1,7),0)) 
  Em=array(0,c(nrow(wide),length(days),8))
  for(i in 1:nrow(wide)) {
    for(j in 1:length(days)){
      state=unlist(wide[i,days[j]])
      Em[i,j,]=emmssionNum[,state]
    }
  }
  return(Em)
}
