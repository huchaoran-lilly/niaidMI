sim_data <-
function(fit,s,days,bin) {
  ret=matrix(0,length(s),length(days))
  colnames(ret)=days
  for(i in 1:length(s)) {
    ret[i,1]=samp(fit$Pri[[ s[i] ]])
    for(t in 2:length(days)) {
      mytran=fit$Tran[[ bin[t-1] ]]
      ret[i,t]=samp(mytran[ret[i,t-1],])
    }
  }
  return(ret)
}
