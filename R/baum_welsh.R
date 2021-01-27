.baum_welsh <-
function(Pri, s, Tran, Em, bin, n_bin=max(bin),
                    n=dim(Em)[1], n_day=dim(Em)[2], n_states=8, n_strata=max(s),
                    tol=1e-6, maxiter=200) {
  
  logLike_new <- 0
  logLike_old <- 1

  nit=0
  while((abs(logLike_new-logLike_old)/abs(logLike_old)) > tol | nit<3) {
    nit=nit+1
    
    if(nit>maxiter) {
      stop("maximum number of iterations reached.")
      # hist(rowSums(log(fb$c)))
      # pos=which(apply(log(fb$c),1, min)< -6.1)
      # 
      # log(fb$c)[pos[2],]
      # cat(logLike_new-logLike_old,  logLike_new,
      #     (abs(logLike_new-logLike_old)/logLike_old),"\n")
      # browser()
    }

    logLike_old <- logLike_new
     
    
    #Forward-Backward application
    fb <- .forward_backward(Pri=Pri, s=s, Tran=Tran, Em=Em, bin=bin,
                           n_day=n_day, n_states=n_states, n_bin=n_bin, n_strata=n_strata)

    
    #Log-Likelihood computation
    fbc=fb$c
    V=fb$beta * Em
    alpha=fb$alpha
    logLike_new <- sum(log(fbc))
    # cat("     ", logLike_new,"\n")
    
    # if(is.na(logLike_new))
    #   browser()
    
    #Update predicted  transitions
    #tran_w = array(0, c(n_bin, n_states,n_states))
    tran_w =replicate(n_bin, matrix(0,n_states,n_states), simplify = FALSE)

    for(i  in 1:n) {
      for(t in 1:(n_day-1)) {
        O1=alpha[i,t,]*Tran[[ bin[t] ]]
        #v1=fb$beta[i,t+1,] * Em[i,t+1,]
        tran_w0=col_mult(O1, V[i,t+1,])
        tran_w[[ bin[t] ]]=tran_w[[ bin[t] ]]+tran_w0/fbc[i,t+1]#sum(tran_w0)
      }
    }
    
    ##Parameters update transition
    
    Tran=lapply(tran_w, function(tmp) {
      norm_const=rowSums(tmp)
      if(any(norm_const==0)) {
        #in this scenario data aparently contains no information 
        #where patient transitions out of a state
        #to avoid divide by zero nan we fix
        #parameter transition out of state should not impact likelihood
        tmp[which(norm_const==0),]=1
        norm_const=rowSums(tmp)
      }
        
      ret=tmp/norm_const
      ret[8,]=c(0,0,0,0, 0,0,0,1)
      return(ret)
    })
    # for(b in 1:n_bin) {
    #   tmp=tran_w[b,,]
    #   tmp=tmp/rowSums(tmp)
    #   tmp[8,]=c(0,0,0,0, 0,0,0,1)
    #   Tran[[b]]=tmp
    # }

    ##Parameters update initial state
    ab1=fb$beta[,1,]*alpha[,1,]
    if(all(s==1)) {
      tmp=colSums(ab1)
      Pri=list(tmp/sum(tmp))
    } else {
      #lapply(split(ab1,s), function(x) colSums(x))
      stop("strata not implemented like this now.") #this should not be possible to happen now
    }

    
  
  }
  
  ## Return and class definition
  bw <- list(Tran=Tran,Pri=Pri,fb=fb, logLike_new=logLike_new, nit=nit)
  return(bw)
}


.forward_backward <-
function(Pri, s, Tran, bin, Em, n=dim(Em)[1], 
                            n_day=dim(Em)[2], n_states=8, n_bin=max(bin),
         n_strata=max(s)) {

  alpha=array(0,c(n, n_day, n_states))
  cm=matrix(0,n, n_day)
  for(i  in 1:n) {
    alpha0 = (Pri[[ s[i] ]] * Em[i,1,])
    cm[i,1] <- sum(alpha0)
    alpha[i,1,] = alpha0 / cm[i,1]
    for(t in 2:n_day) {
      # cat(i,t,"\n")
      alpha0 = (alpha[i,t-1,] %*% Tran[[ bin[t-1] ]]) * Em[i,t,]
      cm[i,t] <- sum(alpha0)
      alpha[i,t, ] <- alpha0 / cm[i,t]
    }
  }
  beta=array(0,c(n, n_day,n_states))
  for(i  in n:1) {
    beta[i,n_day,] =  rep(1,n_states) #/ cm[i,n_day]
    for(t in (n_day-1):1) {
      beta[i,t,] =  (Tran[[bin[t]]] %*% (Em[i,t+1,] * beta[i,t+1,])) / cm[i,t+1]
    }
  }
  return(list(alpha=alpha, beta=beta, c = cm))
}

