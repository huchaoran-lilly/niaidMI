# # bootstrap_param_est_debug <-
# #   function(wide, b, days=paste0("D",1:28), bin=rep(1,length(days)-1),
# #            Em=get_emission(wide, days), tol=1E-6, maxiter=200, silent=FALSE) {
# # 
# # 
# #     if(!is.data.frame(wide))
# #       stop("wide must be a data.frame.")
# # 
# #     if(!is.numeric(bin))
# #       stop("bin must be numeric.")
# #     if(!is.vector(bin))
# #       stop("bin must be a vector")
# #     if(length(bin)!=(length(days)-1))
# #       stop("length(bin) must be the same as length(days)-1.")
# # 
# #     if(!is.numeric(b))
# #       stop("b must be numeric.")
# #     if(!is.vector(b))
# #       stop("b must be a vector")
# #     if(length(b)!=1)
# #       stop("length(b) must be 1.")
# #     if(b<0)
# #       stop("b must be >=0.")
# # 
# #     M = wide[,days]
# #     M = as.matrix(do.call(cbind, lapply(M, function(x) {
# #       x=suppressWarnings(as.numeric(as.character(x)))
# #       # x[x<1 | x>8]=NA
# #       x
# #     })))
# # 
# #     #originally I implemented to stratify baseline. However, not planning to expose this functionality.
# #     #the following line sets up for no strata
# #     s=rep(1, nrow(wide))
# # 
# # 
# #     strt=.get_start(M=M,s=s,bin=bin)
# #     if(!silent) cat("Fiting Model\n")
# #     fit=niaidMI:::.baum_welsh(Pri=strt$Pri, s=s, Tran=strt$Tran, Em=Em, bin=bin, tol=tol,  maxiter=maxiter)
# # 
# #     ret=list(fit=fit, s=s, bin=bin, days=days, Em=Em)
# # 
# #     if(!silent) cat("Fiting Bootstrap Samples:\n")
# #     ret$boot=lapply(1:b, function(i) {
# #       if(!silent) cat("    b=",i)
# #       boot=niaidMI:::.bootstrap_dta(Em=Em, M=M, s=s)
# #       # strt=.get_start(M=boot$M,s=boot$s, bin=bin)
# # 
# #       fit1=niaidMI:::.baum_welsh(Pri=strt$Pri, s=boot$s, Tran=fit$Tran, Em=boot$Em, bin=bin, tol=tol,  maxiter=maxiter)
# #       # fit1=.baum_welsh(Pri=strt$Pri, s=boot$s, Tran=fit$Tran, Em=boot$Em, bin=bin, tol=tol,  maxiter=maxiter)
# #       if(!silent) cat(", iterations=", fit1$nit,"\n")
# #       return(fit1[c("Tran","Pri")])
# #     })
# # 
# #     return(ret)
# #   }
# 
# 
# 
# impute_debug <-
#   function(wide, m, by=NULL, days=paste0("D",1:28), bin=rep(1,length(days)-1), 
#            Em=get_emission(wide, days),
#            listFormatOut=FALSE, tol=1E-6, maxiter=200, silent=FALSE) {
#     
#     
#     if(!is.data.frame(wide))
#       stop("wide must be a data.frame.")
#     
#     if(!is.numeric(bin))
#       stop("bin must be numeric.")
#     if(!is.vector(bin))
#       stop("bin must be a vector")
#     if(length(bin)!=(length(days)-1))
#       stop("length(bin) must be the same as length(days)-1.")
#     
#     if(!is.numeric(m))
#       stop("m must be numeric.")
#     if(!is.vector(m))
#       stop("m must be a vector")
#     if(length(m)!=1)
#       stop("length(m) must be 1.")
#     if(m<2)
#       stop("m must be >=2.")
#     
#     #originally I implemented to stratify baseline. However, not planning to expose this functionality.
#     #the following line sets up for no strata
#     s=rep(1, nrow(wide))
#     
#     
#     if(is.null(by)) { #no by grouping so just go
#       fit=bootstrap_param_est(wide=wide, b=m, bin=bin,
#                               days=days, Em=Em,tol=tol,  
#                               maxiter=maxiter, silent=silent)
#       ### for debug purpose ######################################################
#       #print(fit)
#       ### for debug purpose ######################################################
#       
#       
#       globalcnt<<-0
#       imp=lapply(fit$boot, function(bfit) {
#         globalcnt<<-globalcnt+1
#         .impute_one_debug(Pri=bfit$Pri, s=s, Tran=bfit$Tran, 
#                           Em=Em, 
#                           bin=bin)
#         ### for debug purpose ######################################################
#         
#         ### for debug purpose ######################################################
#         
#       })
#       extravars=wide[,which(!(colnames(wide) %in% days))]
#       
#       ret=mapply(function(impi,i) {
#         colnames(impi)=days
#         data.frame(IMP_ID=i,extravars,impi)
#       }, imp, 1:length(imp), SIMPLIFY = FALSE)
#       
#     } else {
#       
#       if(!is.character(by))
#         stop("'by' must be of type character.")
#       if(!all(by %in% colnames(wide)))
#         stop("All variables listed in by must be in the 'wide' data.frame.")
#       
#       by_str_for_impute=apply(wide[,by, drop=FALSE], 1, function(X) paste(X,collapse="|"))
#       
#       poslst=split(1:nrow(wide),by_str_for_impute)
#       
#       
#       
#       imp=lapply(poslst, function(pos) {
#         # browser()
#         wide_by=wide[pos,]
#         s_by=s[pos]
#         Em_by=Em[pos,,]
#         impute(wide=wide_by, bin=bin, m=m, 
#                days=days, Em=Em_by, 
#                listFormatOut=TRUE, 
#                by=NULL, tol=tol,  
#                maxiter=maxiter,
#                silent=silent)
#       })
#       
#       ret=list()
#       for(i in 1:m)
#         ret[[i]]=do.call(rbind, lapply(imp,function(x) x[[i]] ))
#       
#       ret=lapply(ret, function(x) {
#         rownames(x)=NULL
#         return(x)
#       })
#     }     
#     if(!listFormatOut) {
#       ret=do.call(rbind, ret)
#     }
#     return(ret)
#     
#   }
# 
# 
# 
# 
# .impute_one_debug <-
#   function(Pri, s, Tran, Em, bin, n_bin=max(bin),
#            n=dim(Em)[1], n_day=dim(Em)[2], n_states=8, n_strata=max(s), fix_const=1E-6) {
#     
#     
#     Pri=lapply(Pri,function(x) {
#       x=x+fix_const
#       x/sum(x)})
#     Tran=lapply(Tran,function(x) {
#       x=x+fix_const
#       x/rowSums(x)})
#     
#     fb <- niaidMI:::.forward_backward(Pri=Pri, s=s, Tran=Tran, Em=Em, bin=bin,
#                                       n_day=n_day, n_states=n_states, n_bin=n_bin, n_strata=n_strata)
#     
#     ret=matrix(0, n, n_day)
#     
#     
#     #Log-Likelihood computation
#     fbc=fb$c
#     V=fb$beta * Em
#     alpha=fb$alpha
#     
#     
#     # tran_w =replicate(n_bin, matrix(0,n_states,n_states), simplify = FALSE)
#     for(i  in 1:n) {
#       tmp=fb$beta[i,1,]*alpha[i,1,]
#       p=tmp/sum(tmp)
#       ret[i,1]=niaidMI:::.samp(p)
#       for(t in 1:(n_day-1)) {
#         O1=alpha[i,t,]*Tran[[ bin[t] ]]
#         tran_w0=niaidMI:::col_mult(O1, V[i,t+1,])
#         p=tran_w0[ret[i,t],]
#         p=p/sum(p)
#         ret[i,t+1]=niaidMI:::.samp(p)
#         ### for debug purpose ######################################################
#         #cat(globalcnt, i, t, p, "\n" )
#         ### for debug purpose ######################################################
#       }
#     }
#     
#     return(ret)
#   }