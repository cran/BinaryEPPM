LL.Regression.Binary <-
function(parameter,model.type,model.name,link,ntrials,nsuccess,
               covariates.matrix.p,covariates.matrix.scalef,
               offset.p,offset.scalef,weights,grad.method) {
      nobs <- nrow(covariates.matrix.p) 
      probabilities <- Model.Binary(parameter,model.type,model.name,link,ntrials,
                          covariates.matrix.p,covariates.matrix.scalef,
                          offset.p,offset.scalef)$probabilities 
# Calculation of log likelihood
      vlogl <- rep(0,nobs)
      vlogl <- sapply(1:nobs, function(i) {
         probability <- probabilities[[i]]
         vsuccess    <- nsuccess[[i]] 
         log.prob <- rep(0,length(probability))
         log.prob <- sapply(1:length(probability), function(j) {
                       if ((is.na(probability[j])==TRUE) | (probability[j]<=0)) { 
                                   log.prob[j] <- -1.e+20
                          } else { log.prob[j] <- log(probability[j]) }} ) # end of sapply
         if (is.null(weights)==TRUE) { vlogl[i] <- t(log.prob)%*%vsuccess
                              } else {
            if (is.list(weights)==TRUE) { wk.wts <- weights[[i]]
                                 } else {
               wk.wts <- c(rep(weights[[i]],length(vsuccess)))
                                        } # end of is.list(weights)
         vlogl[i] <- t(wk.wts*log.prob)%*%vsuccess } } ) # end of sapply
      loglikelihood <- sum(vlogl)
   return(loglikelihood) }
