GBprob <-
function(twoparameter,nt) {
      m <- nt + 1
      a <- twoparameter[1] 
      b <- twoparameter[2]
# a = - log(1-p)
      if ((a>0) & ((a<1.e+20) | (a!=Inf))) {
         if (round(b,digits=20)==1) {
            p <- 1 - exp(-a)
            probability <- dbinom(c(0:nt),nt,p,log=FALSE) 
                                          } else {
            vlambda <- c(rep(a,m))
            if (b>0) {            
               vlambda <- vlambda*(c(rep(nt,m)) - c(0:nt))^c(rep(b,m)) } # end if (b>0)
# limiting value for lambda
            lambda.limit <- 745
            vlambda <- sapply(1:m, function(j) 
                 if ((is.finite(vlambda[j])==FALSE) | (vlambda[j]>lambda.limit)) { 
                                        vlambda[j] <- lambda.limit  
                                      } else { vlambda[j] <- vlambda[j] } ) # end of sapply
            probability <- EPPMprob(vlambda) 
                                    } # end of if (b==1)
            } else { 
                 if (a<0) { probability <- c(1,rep(0,nt))
                   } else { probability <- c(rep(0,nt),1) }} # end of if a>0 & a>1.e+20
      return(probability)           }
