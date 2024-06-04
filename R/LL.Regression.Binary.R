LL.Regression.Binary <-
function(parameter,model.type,model.name,link,ntrials,nsuccess,
               covariates.matrix.p,covariates.matrix.scalef,
               offset.p,offset.scalef,weights,grad.method) {
      nobs <- nrow(covariates.matrix.p) 
      probabilities <- Model.Binary(parameter,model.type,model.name,link,ntrials,
                          covariates.matrix.p,covariates.matrix.scalef,
                          offset.p,offset.scalef)$probabilities 

# 1.e-323 = exp(-743.7469) is the smallest probability value not to give -Inf
      small.loglikelihood <- log(1.e-323)
# x = 709 is the largest number not to result in an exp(x) = Inf

# Calculation of log likelihood
      vlogl <- rep(0,nobs)
      for ( i in 1:nobs) {
         probability <- probabilities[[i]]
         vsuccess    <- nsuccess[[i]] 
         log.prob <- ifelse( (((probability>=0) & (probability<1.e-323)) | (is.na(probability)==TRUE)),
                           small.loglikelihood, log(probability))

         if (is.null(weights)==TRUE) { vlogl[i] <- t(log.prob)%*%vsuccess
                              } else {
            if (is.list(weights)==TRUE) { wk.wts <- weights[[i]]
                                 } else {
               wk.wts <- c(rep(weights[[i]],length(vsuccess)))
                                        } # end of is.list(weights)
         vlogl[i] <- t(wk.wts*log.prob)%*%vsuccess } 
                           } # end for loop

      if (is.na(sum(vlogl))==TRUE) { loglikelihood <- nobs*small.loglikelihood
                            } else { loglikelihood <- sum(vlogl) } # end if (is.na(sum(vlogl))==TRUE)
   return(loglikelihood) }
