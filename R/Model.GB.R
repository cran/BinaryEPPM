Model.GB <-
function(parameter,model.name,link,ntrials,
                        covariates.matrix.p,offset.p=c(rep(0,length(ntrials)))) {
#  data as number of trials & number of successes  
   npar <- length(parameter) 
   nobs <- nrow(covariates.matrix.p) 
   vone <- rep(1,nobs) 
   if (model.name=="binomial") { vb <- vone <- rep(1,nobs) 
                        r.parameter <- parameter }
   if (model.name=="EPPM extended binomial") { 
# restricting b to >=0, 0 being Poisson variance, 1 being binomial
                    if (parameter[npar]<0) { parameter[npar] <- 0 }
                    r.parameter <- parameter[1:(npar-1)] 
                    vb <- rep(parameter[npar],nobs) }
   vlp  <- covariates.matrix.p%*%r.parameter + offset.p
# inverse of link function
   vp <- attr(link, which="p")$linkinv(vlp)
   denom <- rep(0,nobs)
   denom <- sapply(1:nobs, function(i) 
             denom[i] <- max(ntrials[[i]]) )
   vmean <- denom*vp
   probabilities <- ntrials 
   if (round(vb[1], digits=20)==1) { va <- - log(vone - vp)
            } else { vonemb <- vone - vb 
                     va     <- (denom^vonemb - (denom - vmean)^vonemb) / vonemb }
   probabilities <- lapply(1:nobs, function(i) 
          probabilities[[i]] <- GBprob(twoparameter=c(va[i],vb[i]),nt=denom[i]) )
   output <- list(model.name=model.name,link=link,parameter=parameter,
                  probabilities=probabilities,
                  Dparameters=data.frame(va,vb))
   return(output) }
