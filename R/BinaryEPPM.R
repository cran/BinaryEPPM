BinaryEPPM <-
function(formula,data,subset=NULL,na.action=NULL,weights=NULL,
                  model.type="p and scale-factor",model.name="generalized binomial",
                  link="cloglog",initial=NULL,method="Nelder-Mead",
                  pseudo.r.squared.type="square of correlation",control=NULL) {

# Checking for correct combinations of model.type and model.name
   if (model.name=="binomial") { 
      if (model.type!="p only") { model.type <- "p only"
         warning("\n","model.type for binomial set to p only") } # end p only
                        } else {
      if ((model.name!="generalized binomial") & (model.name!="beta binomial") &
          (model.name!="correlated binomial")) {
         warning("\n","unknown model.name for this model.type")
         return(object=NULL) } } # end if binomial                       

# Checking method 
   if ((method!="Nelder-Mead") & (method!="BFGS")) {
       warning("\n","unknown function optim method")
       return(object=NULL) }

# Checking that data is data.frame or list.
   if ((is.data.frame(data)==FALSE) & (is.list(data)==FALSE)) { 
      warning("\n","Input data is neither data frame nor list.")
      return(object=NULL) } # end of check data.frame or list

    cl <- match.call()

# link functions
    if (link=="powerlogit") { 
       if (is.null(attr(link, which="power"))==TRUE) {
          attr(link, which="power") <- 1 
                  } else {
          if ((is.finite(attr(link, which="power"))==FALSE) | 
              ((attr(link, which="power")<0)==TRUE)) {
             warning("\n","non finite or <0 power, reset to 1")
             attr(link, which="power") <- 1 }} } # end if link=powerlogit

# Checking for correct link functions
   if ((link!="logit") & (link!="probit") & (link!="cloglog") & 
       (link!="cauchit") & (link!="log") & (link!="loglog") & 
       (link!="doubexp") & (link!="doubrecip") & (link!="powerlogit") & 
       (link!="negcomplog")) {
       warning("\n","unknown link function")
       return(object=NULL) 
           } else {
       if ((link=="doubexp") | (link=="doubrecip") | 
           (link=="powerlogit") | (link=="loglog") |
           (link=="negcomplog")) {
          if (link=="loglog")     { attr(link, which="p") <- loglog() }
          if (link=="doubexp")    { attr(link, which="p") <- doubexp() }
          if (link=="doubrecip")  { attr(link, which="p") <- doubrecip() }
          if (link=="powerlogit") { attr(link, which="p") <- powerlogit(power=attr(link, which="power")) }
          if (link=="negcomplog")  { attr(link, which="p") <- negcomplog() }
                                } else {
          attr(link, which="p") <- make.link(link) } } #end of if link

# as in betareg package handling subset=option 
    if (missing(data)) { data <- environment(formula) }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE 

    FBoth <- Formula(formula) 
    lenFB <- length(FBoth) 
    mf[[1L]] <- as.name("model.frame")

# name of response variable
    wk.name <- attr(FBoth, which="lhs") 
    if (length(wk.name)>1) {
       warning("\n","more than one variable name on lhs of the formula")
       return(object=NULL) }

# data.frame or list

   if (is.data.frame(data)==TRUE) { 

# indicator of whether a data frame (TRUE) or a list (FALSE)
# data frame input
      data.type <- TRUE

      mfBoth    <- model.frame(FBoth,data=data)
      binom.var <- model.part(FBoth,data,lhs=1)
# storing name of numerator in response variable for use 
# with the newdata option in predict.BinaryEPPM
      denom.name <- attr(binom.var, which="names")[2][1]
      resp.var   <- binom.var[[1]]
      n.var      <- binom.var[[2]]
      nvar       <- length(data)
      nobs       <- length(data[[1]])
      if (nvar==1) { covariates  <- NULL
                   } else { covariates <- data } 
      list.data  <- lapply(1:nobs, function(i) {
         if (((resp.var[i]>0)==TRUE) & ((resp.var[i]<n.var[i])==TRUE)) { 
                            c(rep(0,resp.var[i]),1,rep(0,(n.var[i]-resp.var[i])))
                                       } else {
            if ((resp.var[i]==0)==TRUE) { c(1,rep(0,n.var[i])) }
            if ((resp.var[i]==n.var[i])==TRUE) { c(rep(0,n.var[i]),1) } }}) # end of lapply
      p.obs        <- resp.var / n.var
      scalef.obs   <- rep(1,nobs)
      mean.obs     <- resp.var
      variance.obs <- rep(0,nobs)
      vnmax        <- as.vector(n.var)
# storing name numerator of response variable as attribute of vnmax
      attr(vnmax, which="names") <- denom.name
      list.data  <- lapply(1:nobs, function(i) 
                   (resp.var[i]==0)*(c(1,rep(0,n.var[i]))) +
                   (resp.var[i]==n.var[i])*(c(rep(0,n.var[i]),1)) +
                   ((resp.var[i]>0) & (resp.var[i]<n.var[i]))* 
                   (c(rep(0,resp.var[i]),1,rep(0,(n.var[i]-resp.var[i]))))
                         ) # end of lapply

      wkdata <- data.frame(p.obs,scalef.obs,data)

                                   } else { 

# list of frequency distributions input
      data.type <- FALSE

# Checking for name of dependent variable in list of data
      wk.name <- attr(FBoth, which="lhs") 
      nvar <- length(data)
      nobs <- length(data[[1]])
      if (nvar==1) { covariates  <- NULL }

      if ((nvar==1) & (is.list(data[[1]])==TRUE)) { 
         if ((wk.name==names(data)[1])==TRUE) {
              list.set <- TRUE
              list.data <- data[[1]]
                                       } else {          
              warning("\n","single list with list of data is not named ",wk.name)
              return(object=NULL) } # end of if wk.name==names(data)[1]
          } else { list.set <- FALSE
                   for ( i in 1:nvar ) { 
            if (is.list(data[[i]])==TRUE) { 
               if ((wk.name==names(data)[i])==TRUE) { 
                  if (list.set==TRUE)  { 
                     warning("\n","More than one list named ",wk.name," within list of data so",
                             "\n","not clear which is the dependent variable.")
                                         return(object=NULL)
                                } else { list.data <- data[[i]] 
                                         list.set <- TRUE } # end of list.set==TRUE
                                                    } # end of wk.name==names(data)[i]

                                       } else { 
                            if (i==1) { covariates <- data.frame(data[1])
                                        names(covariates[1]) <- names(data[1])                   
                                             } else { 
                               if ((i==2) & ((list.set)==TRUE)) { 
                                               covariates   <- data.frame(data[2])
                                               names(covariates[1]) <- names(data[2]) 
                                                         } else { 
                                               wks <- length(covariates) + 1
                                               covariates <- data.frame(covariates,data[i])
                                               names(covariates[wks]) <- names(data[i])  
                                                           }}} } # end for loop
                                                 } # end if nvar==1 & is.list
      if (list.set==FALSE) { 
         warning("\n","No list named ",wk.name," within list of data.")
         return(object=NULL) } # end of if list.set==FALSE

      mean.obs     <- rep(0,nobs)
      variance.obs <- rep(0,nobs)
      p.obs        <- rep(0,nobs)
      scalef.obs   <- rep(1,nobs)
      vnmax        <- rep(0,nobs)
      for ( i in 1:nobs ) { 
         count        <- list.data[[i]]
         nmax1        <- length(count)
         nmax         <- nmax1 - 1
         vnmax[i]     <- nmax
         cnum         <- 0:nmax
         ncount       <- sum(count)  
         mean.obs[i]  <- t(cnum)%*%count/ncount 
         p.obs[i]     <- mean.obs[i]/nmax
         if (p.obs[i]==0) { p.obs[i] <- 1.e-10 }
         if (p.obs[i]==1) { p.obs[i] <- 1 - 1.e-10 }
         variance.obs[i] <- (t(cnum*cnum)%*%count - ncount*mean.obs[i]*mean.obs[i]) / (ncount - 1)
         scalef.obs[i] <- variance.obs[i] / (mean.obs[i]*(1-p.obs[i]))
         if (is.finite(variance.obs[i])==FALSE) { variance.obs[i] <- 0 } 
         if (is.finite(scalef.obs[i])==FALSE) { scalef.obs[i] <- 1 } 
                    } # end of for loop

# Checking for no covariates 
      if (is.null(covariates)==TRUE) { wkdata <- data.frame(p.obs,scalef.obs) 
                              } else { wkdata <- data.frame(p.obs,scalef.obs,covariates) }

                            }  # end of if is.data.frame

# weights is not null  
   if (is.null(weights)==FALSE) {
      if (is.null(attr(weights, which="normalize"))==TRUE) {
         attr(weights, which="normalize") <- FALSE }

      if (attr(weights, which="normalize")==TRUE) {
         if (is.null(attr(weights, which="norm.to.n"))==TRUE) {
            wkv <- c(rep(0, nobs))
            wkv <- sapply(1:nobs, function(i) { 
                      wkv[i] <- sum(list.data[[i]]) } ) # end of sapply 
            attr(weights, which="norm.to.n") <- sum(wkv) }}

# normalizing weights 

#      if (is.null(attr(weights, which="normalize"))==TRUE) { 
#          attr(weights, which="normalize") <- FALSE }

      if (attr(weights, which="normalize")==TRUE) {
         wkv <- c(rep(0, nobs))
         wkv <- sapply(1:nobs, function(i) { 
                   wkv[i] <- sum(list.data[[i]]) } ) # end of sapply 
         if (is.null(attr(weights, which="norm.to.n"))==TRUE) {
               attr(weights, which="norm.to.n") <- sum(wkv) }
         if (data.type==TRUE) {
            weights <- as.numeric(attr(weights, which="norm.to.n"))*weights/sum(wkv)
                     } else {
            weights <- lapply(1:nobs, function(i) { weights[[i]] <-  
                       as.numeric(attr(weights, which="norm.to.n"))*weights[[i]]/sum(wkv) } ) # end of lapply
                            } } # end if(attr(weights,

				  } # end is.null(weights) 

   mf$data <- wkdata
   mf$resp.var <- p.obs
   mf$formula  <- update(formula(FBoth,lhs=NULL,rhs=1), p.obs ~ . )
# Setting up weights
   if (data.type==FALSE) { mf$weights=vnmax } # end if data.type

# This statement evaluates mf as a data.frame
   wkdata <- eval(mf, parent.frame())

#   evaluating scale-factor model if used 
    if (model.type=="p and scale-factor") { 
          mf.scalef <- mf
          mf.scalef$resp.var <- scalef.obs
       if (lenFB[2]==2) {
          mf.scalef$formula <- update(formula(FBoth,lhs=NULL,rhs=2), scalef.obs ~ . )
                        } else {
          mf.scalef$formula <- update(formula(FBoth,lhs=NULL,rhs=NULL), scalef.obs ~ 1 )
                               } # end if (lenFB[2]==2)

# This statement evaluates mf.scalef as a data.frame
      temp.wkdata <- eval(mf.scalef, parent.frame())
# removing from temp.wkdata the variables common to both data frames, i.e., it & wkdata
      intersection.var <- intersect(names(wkdata), names(temp.wkdata))
      dup.var <- duplicated(c(intersection.var, names(temp.wkdata)))
      wks <- length(intersection.var) + 1
      wke <- length(dup.var)
      dup.var <- dup.var[wks:wke]
      temp.wkdata <- subset(temp.wkdata, select=(dup.var==FALSE))

# merging two data sets i.e., that for mean.obs and that for variance.obs
      wkdata <- merge(wkdata, temp.wkdata, by=0, incomparables=TRUE, sort=FALSE)
                    } # end of if p and scale-factor

# removing the variable (resp.var) and
# removing the variable Row.names added in by the merge operation
# when model type is p and scale-factor 
   wkdata <- subset(wkdata, select=((names(wkdata)!="(resp.var)") &
                                    (names(wkdata)!="Row.names")))

# subsetting mean.obs, variance.obs if subset in operation
   if (is.null(subset)==FALSE) { 
      p.obs        <- p.obs[subset]
      mean.obs     <- mean.obs[subset]
      variance.obs <- variance.obs[subset] 
      scalef.obs   <- scalef.obs[subset] 
      weights      <- weights[subset] 
      vnmax        <- vnmax[subset]
# list.data also needs to be subsetted 
      list.data <- list.data[subset]
# following inserted 7th July 2017 to handle an error resulting when list input is used
      if (data.type==TRUE) { resp.var     <- resp.var[subset] 
         n.var        <- n.var[subset] } } # end of is.null(subset)

# if data frame adding in the two variables resp.var & n.var to the data frame
# for use calculating initial estimates with glm
      if (data.type==TRUE) { wkdata <- data.frame(wkdata, resp.var, n.var) }
      wkdata <- data.frame(wkdata, vnmax)

# Resetting nobs if subset is being used
      if (is.null(subset)==FALSE) { nobs <- length(wkdata[[1]]) }

# Checking if offset.var involved and if it is adding it to workdata so that 
# both variables offset(****) and **** are in the data frame in order to use offsets in glm 
# N.B. DO NOT USE THE NAME offset.*** as an offset name in the data sets input
    offset.var <- NULL
    wks <- length(wkdata)
    for ( i in 1:wks) { nchar <- nchar(names(wkdata[i]))
       if ((nchar>7) & (((substr(names(wkdata[i]), 1, 7)=="offset(")==TRUE) | 
                         (substr(names(wkdata[i]), 1, 7)=="offset.")==TRUE)) { 
             offset.var <- wkdata[i]
             names(offset.var) <- substr(names(wkdata[i]), 8, (nchar-1))
             wkdata <- data.frame(wkdata, offset.var)
                       } } # end of for i

    vone         <- rep(1,nobs)
# link function
    lp.obs <- attr(link, which="p")$linkfun(p.obs)

# Adding lp.obs the linear predictor to wk.data. This is needed for initial estimates
# of parameters when initial not set and links doubexp, doubrecip and power logit used.
    if ((is.null(initial)==TRUE) & ((link=="doubexp") | 
        (link=="doubrecip") | (link=="powerlogit") | (link=="negcomplog"))) {
         wkdata <- data.frame(wkdata, lp.obs) } # end if is.null initial & link

# Changing offsets" names from offset.xxx. to offset(xxx)
    wks <- length(wkdata)
    names(wkdata) <- sapply(1:wks, function(i) {
        nchar <- nchar(names(wkdata[i]))
        if ((nchar>7) & ((substr(names(wkdata[i]), 1, 7)=="offset.")==TRUE)) { 
              new.label <- paste("offset(", substr(names(wkdata[i]), 8, (nchar-1)), ")", sep="") 
              names(wkdata[i]) <- new.label
                    } else { 
              names(wkdata[i]) <- names(wkdata[i]) }
                        } ) # end of sapply

# removing wkdata on exit 
    on.exit(rm(wkdata))

# Checking arguments of function

# Checking for correct model.types
   if ((model.type!="p only") & (model.type!="p and scale-factor")) {
      covariates.matrix.p      <- NULL
      covariates.matrix.scalef <- NULL
      warning("\n","unknown model.type")
      return(object=NULL) 
           } else {
# Checking that formula has 1 lhs and either 1 or 2 rhs
# corresponding to just a formula for p, and a formula for both p
# and scale-factor 
      if (model.type=="p only") {
         if (lenFB[2]==2) { 
            warning("\n","model.type is p only but rhs of formula has two parts to it",
                    "\n","2nd part of rhs of formula is ignored","\n") } # end if lenFB[2]
            FBoth <- update(FBoth,p.obs ~ .) } # end of p only
      if (model.type=="p and scale-factor") {
         if (lenFB[2]==1) { 
# Model for scale-factor set to intercept only.
            FBoth <- update(FBoth,p.obs | scalef.obs ~ . | 1 )
                                     } else {
            FBoth <- update(FBoth,p.obs | scalef.obs ~ . | . ) }
                                       } } # end of if p and scale-factor 

   terms.p <- terms(formula(FBoth,lhs=1,rhs=1))
   if (model.type=="p only") { 
      terms.scalef <- terms(formula( 0 ~ 1 )) 
                             } else { 
      terms.scalef <- terms(formula(FBoth,lhs=2,rhs=2)) } # end of if model.type
   terms.full   <- terms(formula(FBoth))

# Checking that list.data is a list  
if (is.list(list.data)==FALSE) { warning("\n","list.data is not a list") 
   return(object=NULL) }

   nsuccess <- list.data
   ntrials  <- list.data
   ntrials <- lapply(1:nobs, function(i) { nmax1 <- vnmax[i] + 1
                   ntrials[[i]] <- c(rep(vnmax[i],nmax1)) } ) # end of lapply
   names(ntrials) <- NULL

# extracting offsets and model matrices 
   p.mf <- model.frame(formula(FBoth,lhs=1,rhs=1), data=wkdata)

# The following command seems to have problems with handling a variable
# named class. This was the original name of the variable in the Titanic
# data. When changed to pclass there were no problems.
   covariates.matrix.p <- model.matrix(p.mf, data=wkdata) 
   offset.p <- model.offset(p.mf)
   covariates.matrix.scalef <- matrix(c(rep(1,nrow(covariates.matrix.p))),ncol=1)
   offset.scalef <- NULL
   if ((model.type=="p and scale-factor") & (lenFB[2]==2)) {
         scalef.mf <- model.frame(formula(FBoth,lhs=2,rhs=2), data=wkdata)
         covariates.matrix.scalef <- model.matrix(scalef.mf, data=wkdata) 
         offset.scalef <- model.offset(scalef.mf) } # end if p and scale-factor 

   if (is.null(offset.p)==TRUE) { offset.p <- c(rep(0,nobs)) }
   if (is.null(offset.scalef)==TRUE) { offset.scalef <- c(rep(0,nobs)) } 

# Setting up initial estimates of parameters if not input 
   if (is.null(initial)==TRUE) { 
# data.frame input or case data input as list, 
# Using glm to obtain initial estimates 

# setting up new formula in standard form for data.frame input
   if (data.type==TRUE) { 
# changing to standard form of arguments to glm for binomial for data as a data frame
      if (is.null(weights)==TRUE) { 
         FBoth_one <- update(FBoth,cbind(resp.var,(n.var-resp.var)) ~ .) 
         glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=binomial(attr(link, which="p")),
                            data=wkdata)
                           } else { 
         wkdata <- data.frame(wkdata, weights) 
         FBoth_one <- update(FBoth,cbind(resp.var,(n.var-resp.var)) ~ .) 
         glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=binomial(attr(link, which="p")),
                            data=wkdata, weights=weights) }
                     } else {
         weights.p <- c(rep(1,nobs))
         if (is.null(weights)==FALSE) { 
            weights.p <- as.vector(sapply(1:nobs, function(i) { 
                            weights.p[i] <- t(weights[[i]])%*%list.data[[i]] } ) ) } # end of is.null(weights)
         weights.p <- weights.p*vnmax
         wkdata <- data.frame(wkdata, weights.p) 
         FBoth_one <- update(FBoth, p.obs ~ . ) 
         glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=binomial(attr(link, which="p")),
                            data=wkdata, weights=weights.p) 
                     } # end if data.type
      initial.p <- coefficients(glm.Results)

      if (model.type=="p only") {
          if (model.name=="binomial") { parameter <- initial.p
                           names(parameter) <- names(initial.p) }
          if (model.name=="generalized binomial")  { parameter <- c(initial.p,1) 
                           names(parameter) <- c(names(initial.p),"GB parameter") }
          if (model.name=="beta binomial") { parameter <- c(initial.p,0) 
                           names(parameter) <- c(names(initial.p),"beta binomial theta") }
          if (model.name=="correlated binomial") { parameter <- c(initial.p,0) 
                           names(parameter) <- c(names(initial.p),"correlated binomial rho") }
          numpar <- length(parameter)
                                   } # of if model.type p only
      if (model.type=="p and scale-factor") {
            glm.Results <- glm(formula(FBoth,lhs=2,rhs=2),family=gaussian(link="log"),
                                subset=(scalef.obs>0), data=wkdata) 
            initial.scalef <- coefficients(glm.Results)         
            parameter <- c(initial.p,initial.scalef)
            names(parameter) <- c(names(initial.p),names(initial.scalef)) 
            numpar <- length(parameter) } # of if model.type p and scale-factor

                                   } else { # else of initial

# Checking length of input initial against model value
      parameter <- initial
      numpar    <- length(parameter) 
      if (length(initial)!=numpar) { 
         warning("\n","length of initial not equal to number of parameters")
         return(object=NULL) }
      if (is.null(names(initial))==TRUE) { 
         warning("\n","initial has no associated names")
         names(parameter) <- 1:numpar } } # end if is.null(initial)

   start <- parameter

# Checking nobs is the same value as length(list.data)
   wks <- length(ntrials)
   if (wks!=nobs) { 
         warning("\n","number of rows of covariates.matrix.p not equal to length of list.data")
         return(object=NULL) }
   npar.p      <- ncol(covariates.matrix.p)
   npar.scalef <- ncol(covariates.matrix.scalef)

   if (model.type=="p only") {
      npar <- npar.p 
      if (model.name!="binomial") { npar <- npar + 1 }}
   if (model.type=="p and scale-factor") {
      npar <- npar.p + npar.scalef }
# use of the waldtest function with initial estimates used for
# the original call to BinaryEPPM causes this error if 
# (is.null(initial)==TRUE) is not included 
   if ((numpar!=npar) & (is.null(initial)==TRUE)) { 
      warning("\n","number of parameters error") 
      return(object=NULL) }

# link function for mean
      r.parameter.p <- rep(0,npar.p) 
      r.parameter.p <- parameter[1:npar.p] 
      lp.p <- covariates.matrix.p%*%r.parameter.p + offset.p
# link function for variance 
      if (model.type=="p and scale-factor") { 
         r.parameter.scalef <- rep(0,npar.scalef) 
         wks <- npar.p + 1
         r.parameter.scalef <- parameter[wks:npar] 
         lp.scalef <- covariates.matrix.scalef%*%r.parameter.scalef + offset.scalef } # end of if statement

      if ((pseudo.r.squared.type!="square of correlation") & (pseudo.r.squared.type!="R squared") & 
          (pseudo.r.squared.type!="max-rescaled R squared")) {
            warning("\n","unknown argument for pseudo.r.squared.type") 
            return(object=NULL) } # end of if pseudo.r.squared.type

# Checking grad.method if method is BFGS
      if (method=="BFGS") { 
         if (is.null(attr(method,which="grad.method"))==TRUE) { 
               attr(method,which="grad.method") <- "simple" }
         if ((attr(method,which="grad.method")!="simple") & 
             (attr(method,which="grad.method")!="Richardson")) {
            warning("\n","unknown gradient method with method BFGS so reset to simple")
               attr(method,which="grad.method") <- "simple" }
               grad.method <- attr(method,which="grad.method")
                           } else {
               grad.method <- NULL
            } # end if method=BFGS

# Checking for initial estimates being in legitimate range       
      if (is.null(initial)==FALSE) {
         initial.loglikelihood <- LL.Regression.Binary(parameter,model.type,model.name,link,
                                     ntrials,nsuccess,covariates.matrix.p,covariates.matrix.scalef,
                                     offset.p,offset.scalef,weights,grad.method)
         if (initial.loglikelihood<=-1.e+20) { 
            warning("\n","initial estimates give a log likelihood outside legitimate range")
            return(object=NULL) } } # end if is.null(initial)

# Setting defaults and checking control parameters for optim
   if (is.null(control)==TRUE) { 
      control=list(fnscale=-1,trace=0,maxit=1000,abstol=1e-8,reltol=1e-8,
                   alpha=1.0,beta=0.5,gamma=2.0,REPORT=10) 
                               } else { 
# Control parameters related to Nelder-Mead except for REPORT which
# is for BFGS. Those related only to other methods of optim are considered to 
# be unknown and ignored.
     ncontrol <- length(control)
     wk.control <- list(fnscale=-1,trace=0,maxit=1000,abstol=1e-8,reltol=1e-8,
                        alpha=1.0,beta=0.5,gamma=2.0,REPORT=10) 
     for ( i in 1:ncontrol) { 
         if ((names(control[i])!="fnscale") & (names(control[i])!="trace") & 
             (names(control[i])!="maxit") & (names(control[i])!="abstol") & 
             (names(control[i])!="reltol") & (names(control[i])!="alpha") & 
             (names(control[i])!="beta") & (names(control[i])!="gamma") & 
             (names(control[i])!="REPORT")) { 
            wkname <- names(control[i])
            warning("\n","Argument in control is unknown so ignored.")
                                  } # end of if
         if (names(control[i])=="fnscale") { wk.control$fnscale <- control[[i]] }
         if (names(control[i])=="trace")   { wk.control$trace   <- control[[i]] }
         if (names(control[i])=="maxit")   { wk.control$maxit   <- control[[i]] }
         if (names(control[i])=="abstol")  { wk.control$abstol  <- control[[i]] }
         if (names(control[i])=="reltol")  { wk.control$reltol  <- control[[i]] }
         if (names(control[i])=="alpha")   { wk.control$alpha   <- control[[i]] }
         if (names(control[i])=="beta")    { wk.control$beta    <- control[[i]] }
         if (names(control[i])=="gamma")   { wk.control$gamma   <- control[[i]] }
         if (names(control[i])=="REPORT")  { wk.control$REPORT  <- control[[i]] }
                             } # end of for loop
     control <- wk.control
                                   } # end of is.null(control)

# Fitting model using optim, options Nelder-Mead with gr=NULL, BFGS with gr=gradient function

       if (length(parameter)==1) {
# Using estimate of p from the binomial model
          converged <- TRUE 
          attr(converged, which="code") <- NA
          intercept.loglikelihood <- LL.Regression.Binary(parameter=coefficients(glm.Results),
                                           model.type=model.type,model.name=model.name,link=link,
                                           ntrials=ntrials,nsuccess=list.data,
                                           covariates.matrix.p=matrix(c(rep(1,nobs)),ncol=1),
                                           covariates.matrix.scalef=NULL,
                                           offset.p=offset.p,offset.scalef=NULL,
                                           weights=weights,grad.method=NULL) 
          wk.optim <- list(par=coefficients(glm.Results), value=intercept.loglikelihood, 
                           counts=NA, convergence=0, message=NULL)   
                     } else {
          if (method=="Nelder-Mead") { gr.fun <- NULL 
                                 } else { gr.fun <- LL.gradient }              
          wk.optim <- optim(parameter,fn=LL.Regression.Binary,gr=gr.fun,
                            model.type,model.name,link,ntrials,nsuccess,
                            covariates.matrix.p,covariates.matrix.scalef,
                            offset.p,offset.scalef,weights,grad.method,
                            method=method,control=control,hessian=FALSE) 
          if (wk.optim$convergence==0) { converged <- TRUE 
                               } else { converged <- FALSE } 
          attr(converged, which="code") <- wk.optim$convergence 

# The parameter estimates returned from a new call to Model.Binary being different
# from those used as input to the new call to Model.Binary is an indication
# that a limit has been reached. In Model.Binary the distribution parameters are 
# set to the limiting values when those limiting values are reached. 
# Checking for whether scale factor pdistribution parameters are on their limits
# using a comparison of parameter values.
          iprt <- 0
          for ( i in 1:100) { check.optim.par <- 
             Model.Binary(wk.optim$par,model.type,model.name,link,ntrials,
                          covariates.matrix.p,covariates.matrix.scalef,
                          offset.p,offset.scalef)$parameter
             if (max(abs(check.optim.par-wk.optim$par))>1.e-12) { iprt <- 1
                wk.optim$par <- check.optim.par 
                                        } else { break } } # end of for loop
          wk.optim$value <- LL.Regression.Binary(parameter=wk.optim$par,
                               model.type,model.name,link,ntrials,nsuccess,
                               covariates.matrix.p,covariates.matrix.scalef,
                               offset.p,offset.scalef,weights,grad.method)
          if (iprt>0) {
             if (model.name=="generalized binomial") { 
                if (model.type=="p only") { warning("The parameter b=0 showing that the variance", 
                                    "\n","has reached the Poisson boundary hence its se is set to NA.")
                         } else { warning("The parameter b=0 for some observations showing that", 
                                          "\n","the variance has reached the Poisson boundary.") }
                                } # end of if model
             if (model.name=="beta binomial") { 
                if (model.type=="p only") { warning("The value of theta is less than the lower limit", 
                                                    "\n","for some observations hence its se is set to NA.") 
                         } else { warning("The values of theta for some observations are less","\n", 
                                          "than the lower limits.","\n") } } # end of if model
             if (model.name=="correlated binomial") { 
                if (model.type=="p only") { warning("The value of rho is outside the lower to upper limit", 
                                                    "\n","range for some observations hence its se is set to NA.") 
                         } else { warning("The values of rho for some observations are",
                                          "\n","outside the lower to upper limit range.") } } } # end of if iprt==1
                                 } # end of if length(parameter)=1 
          names(wk.optim$par) <- names(parameter) 

          nobs <- nrow(covariates.matrix.p) 
          p.par        <- rep(0,nobs)
          scalef.par   <- rep(1,nobs)
          scalef.limit <- rep(0,nobs)
# Calculation of p and means from parameter estimates and design matrices
          npar.p      <- ncol(covariates.matrix.p)
          r.parameter.p <- rep(0,npar.p) 
          r.parameter.p <- wk.optim$par[1:npar.p] 

# inverse of link function
          p.par <- attr(link, which="p")$linkinv(lp.p)
          if (model.type=="p only") { 
             npar  <- npar.p + 1
             if (model.name=="binomial") { 
                wkv.coefficients <- list(p.est=wk.optim$par[1:npar.p], scalef.est=NULL)
                                   } else { 
                wkv.coefficients <- list(p.est=wk.optim$par[1:npar.p], scalef.est=wk.optim$par[npar]) }
                                    } # if model.type

          if (model.type=="p and scale-factor") { 
# Calculation of scale-factors and variances from parameter estimates and design matrices
             npar.scalef <- ncol(covariates.matrix.scalef)
             npar <- npar.p + npar.scalef
             wks <- npar.p + 1
             wkv.coefficients <- list(p.est=wk.optim$par[1:npar.p], scalef.est=wk.optim$par[wks:npar]) 
# inverse link for scale-factor 
             r.parameter.scalef <- rep(0,npar.scalef) 
             denom <- rep(0,nobs)
             denom <- sapply(1:nobs, function(i) denom[i] <- max(ntrials[[i]]) )
# modeling scalefactor
             vsf  <- as.vector(covariates.matrix.scalef%*%r.parameter.scalef + offset.scalef)
             scalef.par <- exp(vsf)
                                              } # if model.type

# model.hessian from hessian from numDeriv method Richardson
          model.hessian <- hessian(LL.Regression.Binary,
                x=wk.optim$par,method="Richardson",method.args=list(r=6,d=0.001), 
                model.type=model.type,model.name=model.name,link=link,
                ntrials=ntrials,nsuccess=nsuccess,
                covariates.matrix.p=covariates.matrix.p,
                covariates.matrix.scalef=covariates.matrix.scalef,
                offset.p=offset.p,offset.scalef=offset.scalef,weights=weights,
                grad.method=grad.method)

# checking condition number of vcov and inverting plus square rooting to get variance/covariance matrix
          deter <- det(model.hessian)
          wk.npar <- nrow(model.hessian)
          if ((is.finite(deter)==FALSE) | (deter==0)) { vcov <- matrix(c(rep(NA,(wk.npar*wk.npar))),ncol=wk.npar)
                        } else { if (wk.npar==1) { vcov <- -1/model.hessian
                                          } else { condition <- rcond(model.hessian)
# function rcond is from the package Matrix and gives the (reciprocal) condition number
# near 0 is ill-conditioned, near 1 is well conditioned
                                    if (condition>1e-16) {
                                        vcov <- - solve(model.hessian)
                                           } else { vcov <- matrix(c(rep(NA,(wk.npar*wk.npar))),ncol=wk.npar) }}} 
# Setting up row and column names for vcov
          colnames(vcov) <- rownames(vcov) <- names(wk.optim$par[1:wk.npar]) 
          probabilities <- Model.Binary(wk.optim$par,model.type,model.name,link,ntrials,
                                 covariates.matrix.p,covariates.matrix.scalef,
                                 offset.p,offset.scalef)$probabilities 

# Calculate the fitted value and residuals vectors using
# means and variances from predicted probabilities
          mean.prob     <- rep(0,nobs)
          variance.prob <- rep(0,nobs)
          p.prob        <- rep(0,nobs)
          scalef.prob   <- rep(0,nobs)
          for ( i in 1:nobs) { probability      <- probabilities[[i]]
                               nmax             <- vnmax[i] 
                               vid              <- c(0:nmax)
                               fmean            <- t(probability)%*%vid 
                               mean.prob[i]     <- fmean 
                               variance.prob[i] <- t(probability)%*%((vid-as.vector(fmean))^2) 
                               p.prob[i]        <- mean.prob[i]/nmax
                               scalef.prob[i]   <- variance.prob[i] / (mean.prob[i]*(1-p.prob[i]))
                             } # end of for loop

# testing to see if logistic regression in which case the default pseudo R-squared must be either
# R squared or max-rescaled R squared
          wks <- sum(vnmax)
          if ((sum(vnmax)==length(vnmax)) & (pseudo.r.squared.type=="square of correlation")) {       
             warning("logistic regression but pseudo.r.squared.type is square of correlation")
                  } else {
# calculation of pseudo R squared
          if (nobs>1) {
             if (pseudo.r.squared.type=="square of correlation") {
# as in betareg
                if (npar.p>0) {
                   eta <- as.vector(covariates.matrix.p %*% wkv.coefficients$p.est + offset.p)
                        } else {
                   eta <- offset.p }
# link function
                lp.obs <- attr(link, which="p")$linkfun(p.obs)
                eta <- sapply(1:nobs, function(i) {
                    if (abs(eta[i])==Inf) { eta[i] <- NA 
                            } else { eta[i] <- eta[i] }} )
                lp.obs <- sapply(1:nobs, function(i) {
                    if (abs(lp.obs[i])==Inf) { lp.obs[i] <- NA
                            } else { lp.obs[i] <- lp.obs[i] }} )
                if ((sd(eta, na.rm=TRUE)>0) & (sd(lp.obs, na.rm=TRUE)>0)) { 
                     pseudo.r.squared <- cor(eta, lp.obs, use="complete.obs")^2 
                              } else {
#                     warning("One or both of observed or fitted linear predictor has 0 sd.")
                     pseudo.r.squared <- NA } } # end of if square of correlation
            if ((pseudo.r.squared.type=="R squared") | (pseudo.r.squared.type=="max-rescaled R squared")) {
# calculation of pseudo R squared as the generalized coefficient of determination 
# of Cox and Snell (1989) and Nagelkerke (1991).

# obtain log likelihood of intercept only models for a binomial
# Using glm to obtain estimate of p for a binomial model with intercept only

# setting up new formula in standard form for data.frame input
# changing to standard form of arguments to glm for binomial for data as a data frame
            if (data.type==TRUE) { 
               if (is.null(weights)==TRUE) { 
                  FBoth_two <- update(FBoth,cbind(resp.var,(n.var-resp.var)) ~ 1 ) 
                  glm.Results <- glm(formula(FBoth_two,lhs=1,rhs=1),family=binomial(attr(link, which="p")),
                                     data=wkdata) 
                                    } else {
                  wkdata <- data.frame(wkdata, weights) 
                  FBoth_two <- update(FBoth,cbind(resp.var,(n.var-resp.var)) ~ 1 ) 
                  glm.Results <- glm(formula(FBoth_two,lhs=1,rhs=1),family=binomial(attr(link, which="p")),
                                     data=wkdata, weights=weights) } # end of is.null(weights)
                             } else {
               FBoth_two <- update(FBoth, p.obs ~ 1 ) 
               glm.Results <- glm(formula(FBoth_two,lhs=1,rhs=1),family=binomial(link=attr(link, which="p")),
                                     data=wkdata, weights=weights.p) 
                             } # end if data.type 

            intercept.loglikelihood <- LL.Regression.Binary(parameter=coefficients(glm.Results),
                                         model.type="p only",model.name="binomial",link=link,
                                         ntrials=ntrials,nsuccess=list.data,
                                         covariates.matrix.p=matrix(c(rep(1,nobs)),ncol=1),
                                         covariates.matrix.scalef=NULL,
                                         offset.p=offset.p,offset.scalef=offset.scalef,
                                         weights=weights,grad.method=NULL) 
            if (data.type==TRUE) { wks <- sum(vnmax)
                          } else { wkv <- c(rep(0,nobs))
                                   wkv <- sapply(1:nobs, function(i)
                                          wkv[i] <- sum(list.data[[i]]) )                                  
                                   wks <- sum(wkv*vnmax) } # end if data.type
            pseudo.r.squared <- (1-exp(2*(intercept.loglikelihood-wk.optim$value)/wks))
            if (pseudo.r.squared.type=="max-rescaled R squared") {
                pseudo.r.squared <- pseudo.r.squared / (1-exp(2*(intercept.loglikelihood)/wks)) }
                                                      } # end of if R square or max-rescaled
                   } else { pseudo.r.squared <- NA } } # end of if ((sum(vnmax)==length(vnmax))
         attr(pseudo.r.squared, which="names") <- pseudo.r.squared.type

# For the list data type (frequency distribution of data) the degrees of freedom
# for null and residual need to be set to those of the data frame data type

        if (data.type==TRUE) { total.ninlist <- nobs
                      } else {
           vninlist <- c(rep(0,length(list.data)))
           vninlist <- sapply(1:length(list.data), function(ilist) 
                             vninlist[ilist] <- sum(list.data[[ilist]]) )
           total.ninlist <- sum(vninlist)
                             } # end of is.data.frame

         object <- list(data.type=data.type, list.data=list.data, call=cl, 
                       formula=formula, model.type=model.type, model.name=model.name, 
                       link=link, covariates.matrix.p=covariates.matrix.p,
                       covariates.matrix.scalef=covariates.matrix.scalef,
                       offset.p=offset.p,offset.scalef=offset.scalef,
                       coefficients=wkv.coefficients,loglik=wk.optim$value,vcov=vcov,
                       n=nobs, nobs=nobs, df.null=total.ninlist, df.residual=(total.ninlist-length(wk.optim$par)),
                       vnmax=vnmax, weights=weights,converged=converged, method=method,
                       pseudo.r.squared=pseudo.r.squared, start=start, optim=wk.optim,
                       control=control, fitted.values=p.prob, y=p.obs,
                       terms=list(p=terms.p,scale.factor=terms.scalef,full=terms.full)) 

     attr(object, "class") <- c("BinaryEPPM")

     return(object) }
