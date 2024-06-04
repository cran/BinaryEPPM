Model.JMVGB <-
function(parameter,model.name,link,ntrials,
                        covariates.matrix.p,covariates.matrix.scalef,
                        offset.p=c(rep(0,length(ntrials))),
                        offset.scalef=c(rep(0,length(ntrials)))) {
#  data as number of trials & number of successes  
   numpar <- length(parameter) 
   nobs   <- nrow(covariates.matrix.p)
   npar.p      <- ncol(covariates.matrix.p)
   npar.scalef <- ncol(covariates.matrix.scalef)
   npar <- npar.p + npar.scalef
   if (npar!=numpar) {
      warning("\n","no. of parameters not equal to sum of no. of columns mean & variance matrices") }
   va <- rep(0,nobs) 
   vb <- rep(0,nobs) 
   vone <- rep(1,nobs) 
#  model can only be EPPM extended binomial
# modeling p
   r.parameter <- rep(0,npar.p) 
   r.parameter <- parameter[1:npar.p]
   vlp <- covariates.matrix.p%*%r.parameter + offset.p
   vone <- rep(1,nobs) 
   vp <- attr(link, which="p")$linkinv(vlp)
   denom <- rep(0,nobs)
   denom <- sapply(1:nobs, function(i) 
             denom[i] <- max(ntrials[[i]]) )
   vmean     <- denom*vp
# modeling scalefactor
   nparm1      <- npar.p + 1
   r.parameter <- rep(0,npar.scalef) 
   r.parameter <- parameter[nparm1:npar]
   vsf  <- as.vector(covariates.matrix.scalef%*%r.parameter + offset.scalef)
   vscalefact <- exp(vsf)
   probabilities <- ntrials 
   vlower <- (vone/(vone-vp)-vscalefact)
   vlower <- sapply(1:nobs, function(i) {
               if (is.na(vlower[i])==TRUE) { vlower[i] <- 1.e+20
                     } else {  
                  if (vlower[i]<=0) { vlower[i] <- 1.e-10
                        } else { vlower[i] <- vlower[i] }}} ) # end of sapply
   vupper <- -vscalefact
   for ( i in 1:nobs) { 
      uniroot.output <- uniroot(function(b,p,scalef) { 
                     if ((p<=0) | (p>=1)) {
                        if (p<=0) { fvalue <- 1 
                                  } else { fvalue <- (1+scalef)/(2*scalef) }
                                          } else { 
                        wk.2bm1 <- 2*b - 1
                        fvalue <- ((1-p)^wk.2bm1 - 1)/(-wk.2bm1*p) - scalef }
                        return(fvalue) },p=vp[i],scalef=vscalefact[i],
                        lower=0,upper=1.e+10,f.lower=vlower[i],f.upper=vupper[i],
                        extendInt="no",check.conv=TRUE,tol=1.e-10,maxiter=100,
                        trace=TRUE)
      vb[i] <- uniroot.output$root

# restricting maximum variance to Poisson i.e., variance=mean, b=0
#      if (vb[i]<0) { vb[i] <- 0 }
# There is no need to check that b=0 as the lower limit for it in uniroot is 0
      if (round(vb[i],digits=20)==1) { 
                    probabilities[[i]] <- dbinom(c(0:denom[i]),denom[i],vp[i],log=FALSE)
           } else { onemb <- 1 - vb[i]
                    va[i] <- (denom[i]^onemb - (denom[i] - vmean[i])^onemb) / onemb
                    probabilities[[i]] <- GBprob(twoparameter=c(va[i],vb[i]),nt=denom[i]) }
                     } # end of for loop

# returning estimates of vb to estimates of the parameters of the 
# linear predictor for scalefact
   v2bm1 <- c(rep(2,nobs))*vb - vone
   wk.vscalefact <- ((vone-vp)**v2bm1 - vone) / (-v2bm1*vp)
# scale-factor can be less than 0, so setting scale-factor = 1.e-10 when that is so
   wk.vscalefact <- ifelse( ((is.na(wk.vscalefact)==TRUE) | (wk.vscalefact<=0)), 1.e-10, wk.vscalefact)
   wk.r.parameter <- qr.solve(covariates.matrix.scalef, (log(wk.vscalefact) - offset.scalef))
   parameter[nparm1:npar] <- wk.r.parameter 

   output <- list(model.name=model.name,link=link,parameter=parameter,
                  probabilities=probabilities,
                  Dparameters=data.frame(va,vb))
   return(output) }
