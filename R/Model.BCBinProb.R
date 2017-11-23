Model.BCBinProb <-
function(parameter,model.type,model.name,link,ntrials,
                        covariates.matrix.p,covariates.matrix.scalef= 
                            matrix(c(rep(1,nrow(covariates.matrix.p))),ncol=1),
                        offset.p=c(rep(0,length(ntrials))),
                        offset.scalef=c(rep(0,length(ntrials)))) {
#  data as number of trials & number of successes  
   twoparameter <- rep(0,2)  
#  the first element of twoparameter is p
#  the second element of twoparameter is theta or rho 
   numpar <- length(parameter) 
   nobs   <- nrow(covariates.matrix.p)
   npar.p      <- ncol(covariates.matrix.p)
   npar.scalef <- ncol(covariates.matrix.scalef)
   npar <- npar.p + npar.scalef
   if (numpar!=npar) {
      cat("\n","no. of parameters not equal to sum of no of columns mean & variance matrices","\n") }

# modeling p 
   r.parameter <- rep(0,npar.p) 
   r.parameter <- parameter[1:npar.p]
   if (npar.p==0) { vlp <- offset.p
           } else { 
      vlp <- as.vector(covariates.matrix.p%*%r.parameter + offset.p) }
   vone <- rep(1,nobs) 
   vp <- attr(link, which="p")$linkinv(vlp)
   denom <- rep(0,nobs)
   denom <- sapply(1:nobs, function(i) 
             denom[i] <- max(ntrials[[i]]) )
   vmean <- denom*vp
   probabilities <- ntrials 
# modeling scalefactor 
   if (model.type=="p and scale-factor") {
      nparm1      <- npar.p + 1
      r.parameter <- rep(0,npar.scalef) 
      r.parameter <- parameter[nparm1:npar]
      if (npar.scalef==0) { vsf <- offset.scalef
           } else { 
         vsf  <- as.vector(covariates.matrix.scalef%*%r.parameter + offset.scalef) }
      vscalefact <- exp(vsf) } # end if model.type
# calculating limits beta and correlated binomial
   lower.limit <- rep(0,nobs) 
   if (model.name=="beta binomial") { 
      lower.limit <- sapply(1:nobs, function(i) {
         if ((vp[i]>0) & ((1-vp[i])>0)) {
# limits as in Smith & Ridout A Remark on Algorithm AS 189 
            if (denom[i]==1)   { lower.limit[i] <- -1 
                               } else {
               if (vp[i]>0.5) { lower.limit[i] <- (vp[i] - 1)/(denom[i]-1) 
                       } else { lower.limit[i] <- -vp[i]/(denom[i]-1) }}
                                        } else {
            lower.limit[i] <- 0 }} ) # end of sapply
                                               } # end of if beta binomial

# calculating limits for correlated binomial
   if (model.name=="correlated binomial") { 
      lower.limit <- sapply(1:nobs, function(i) {
         if ((vp[i]>0) & ((1-vp[i])>0)) {
# limits on rho given in Kupper & Haseman Biometrics (1978) 
            if (denom[i]==1)   { lower.limit[i] <- -1 
                               } else {
               if (vp[i]>0.5) { lower.limit[i] <- - 2*(1-vp[i])/(denom[i]*(denom[i]-1)*vp[i]) 
                       } else { lower.limit[i] <- - 2*vp[i]/(denom[i]*(denom[i]-1)*(1-vp[i])) }}
                                        } else {
            lower.limit[i] <- 0 }} ) # end of sapply
      upper.limit <- rep(0,nobs)
      upper.limit <- sapply(1:nobs, function(i) {
         if ((vp[i]>0) & ((1-vp[i])>0)) {
# limits on rho given in Kupper & Haseman Biometrics (1978) 
            if (denom[i]==1) { upper.limit[i] <- 1 
                      } else {
                gamma0 <- min((c(0:denom[i]) - (denom[i] - 1)*vp[i] - 0.5)**2) 
                upper.limit[i] <- 2*vp[i]*(1-vp[i])/((denom[i]-1)*vp[i]*(1-vp[i]) + 0.25 - gamma0) }
                                        } else {
            upper.limit[i] <- 0 }} ) # end of sapply
                                          } # end of if correlated binomial

   if (model.type=="p only") { 

      vtheta <- rep(parameter[npar],nobs)
      if (model.name=="beta binomial") { 
         if (parameter[npar]<max(lower.limit)) { 
            parameter[npar] <- max(lower.limit)
            vtheta <- rep(max(lower.limit),nobs) }} # end of if beta binomial

      if (model.name=="correlated binomial") { wks <- parameter[npar]
         if (parameter[npar]<max(lower.limit)) { wks <- max(lower.limit) }
         if (parameter[npar]>min(upper.limit)) { wks <- min(upper.limit) }
         parameter[npar] <- wks 
         vtheta <- rep(wks,nobs) } # end of if correlated binomial
                             } # end of if model.type p only

   if (model.type=="p and scale-factor") {

# scale-factor can be undetermined when p=1
# setting scale-factor = 1 making theta = 0
      vscalefact <- sapply(1:nobs, function(i) {
          if (is.finite(vscalefact[i])==FALSE) { vscalefact[i] <- 1
              } else { vscalefact[i] <- vscalefact[i] } } ) # end of sapply
# scale-factor can be less than 0, so setting scale-factor = 1.e-10 when that is so
      vscalefact <- sapply(1:nobs, function(i) {
          if (vscalefact[i]<=0) { vscalefact[i] <- 1.e-10
              } else { vscalefact[i] <- vscalefact[i] } } ) # end of sapply
      vvariance <- vmean*(vone-vp)*vscalefact
      vtheta <- rep(0, nobs)
      vtheta <- sapply(1:nobs, function(i) {
           if (denom[i]==1) { vtheta[i] <- vtheta[i]
                            } else { 
           vtheta[i] <- (vscalefact[i] - 1)/(denom[i]-1) } } ) # end of sapply

# converting rho to theta for beta binomial so for the
# correlated binomial theta is actually rho
      if (model.name=="beta binomial") { 
         vtheta <- vtheta/(vone-vtheta) 
         vtheta <- sapply(1:nobs, function(i) {
             if (vtheta[i]<lower.limit[i]) { vtheta[i] <- lower.limit[i]
                 } else { vtheta[i] <- vtheta[i] } } ) # end of sapply
         wk.vscalefact <- vtheta*(denom-vone)/(vone+vtheta) + vone
                                       } # end of if model.name beta binomial
      if (model.name=="correlated binomial") {
         vtheta <- sapply(1:nobs, function(i) {
             if (vtheta[i]<lower.limit[i]) { vtheta[i] <- lower.limit[i]
                 } else { vtheta[i] <- vtheta[i] } } ) # end of sapply
         vtheta <- sapply(1:nobs, function(i) {
             if (vtheta[i]>upper.limit[i]) { vtheta[i] <- upper.limit[i]
                 } else { vtheta[i] <- vtheta[i] } } ) # end of sapply
         wk.vscalefact <- vtheta*(denom-vone) + vone
                                       } # end of if model.name correlated binomial
# scale-factor can be less than 0, so setting scale-factor = 1.e-10 when that is so
      wk.vscalefact <- sapply(1:nobs, function(i) {
          if (wk.vscalefact[i]<=0) { wk.vscalefact[i] <- 1.e-10
                            } else { wk.vscalefact[i] <- wk.vscalefact[i] } } ) # end of sapply

# returning estimates of the scale factor to estimates of the parameters 
# of the linear predictor        
      wk.r.parameter <- qr.solve(covariates.matrix.scalef, (log(wk.vscalefact) - offset.scalef))
      parameter[nparm1:npar] <- wk.r.parameter 
                                       } # end if model.type p and scale-factor
      vtheta <- as.vector(vtheta)
      names(vtheta) <- NULL

#  Calculating non factorial part of probabilities
   probabilities <- lapply(1:nobs, function(i) { nt <- denom[i] 
      twoparameter[1] <- vp[i] 
      twoparameter[2] <- vtheta[i]
      if (model.type=="p and scale-factor") { twoparameter[2] <- vtheta[i] } 
      if (model.name=="beta binomial") { 
         wk.output <- BBprob(twoparameter,nt) } # end of beta binomial
      if (model.name=="correlated binomial") { 
         wk.output  <- CBprob(twoparameter,nt) } # end of correlated binomial
      probabilities[[i]] <- wk.output$probability } ) # end of lapply

#  Calculating factorials and applying to probabilities
#  but only for p not equal to 0 or 1
   probabilities <- lapply(1:nobs, function(i) { 
      if ((vp[i]>0) & ((1-vp[i])>0)) {
         nt <- denom[i] 
         m <- nt + 1
         vfac <- choose(c(rep(nt,m)),c(0:nt))
         probabilities[[i]] <- probabilities[[i]]*vfac
                                     } else {
         probabilities[[i]] <- probabilities[[i]] }} ) # end of lapply
#  For the correlated binomial probabilities < 0 can sometimes be produced.
#  These are relatively small i.e., of the order of 10^-15 or smaller, so rounding 
#  such probabilities to 0.
   probabilities <- lapply(1:nobs, function(i) { 
      probabilities[[i]] <- sapply(1:length(probabilities[[i]]), function(j) { 
         if (probabilities[[i]][j]<0) { 
            probabilities[[i]][j] <- 0
                                         } else {
            probabilities[[i]][j] <- probabilities[[i]][j] } } ) # end of sapply
                                         } ) # end of lapply
#  Total probability can be slightly different from 1 therefore adjusting to total
   probabilities <- lapply(1:nobs, function(i) { 
      if ((vp[i]>0) & ((1-vp[i])>0)) {
         total.prob <- sum(as.vector(probabilities[[i]]))
# rounding errors can cause the total probability to be slightly different from 1
         probabilities[[i]] <- probabilities[[i]] / total.prob
                                      } else {
         probabilities[[i]] <- probabilities[[i]] } } ) # end of lapply

      if (model.type=="p only") { vtheta <- rep(max(vtheta),nobs) }
      if (model.name=="beta binomial") { 
          output <- list(model.name=model.name,link=link,parameter=parameter,
                         probabilities=probabilities,
                         Dparameters=data.frame(vp,vtheta,lower.limit)) }
      if (model.name=="correlated binomial") { 
          output <- list(model.name=model.name,link=link,parameter=parameter,
                         probabilities=probabilities,
                         Dparameters=data.frame(vp,vrho=vtheta,lower.limit,
                                     upper.limit)) }
   return(output) }
