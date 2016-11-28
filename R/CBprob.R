CBprob <-
function(twoparameter,nt) {
   p   <- twoparameter[1] 
   q   <- 1 - p 
   rho <- twoparameter[2]
#  m should equal nt + 1 
   m  <- nt + 1
   probability <- rep(0,m)

   if ((p>0) & (q>0)) {
        logp <- log(p)
        logq <- log(q)
        np  <- nt*p 
        np2 <- np*p 
        twopm1  <- 2*p - 1 
        twop2q2 <- 2*(p*q)**2 
        probability <- sapply(1:m, function(i) { 
           im1 <- i - 1 
           wks <- 1 + rho*( (im1-np)**2 + im1*twopm1 - np2 )/twop2q2 
           if (wks>=0) { 
              probability[i] <- wks*exp(im1*logp + (nt-im1)*logq)
                       } else { probability[i] <- 0 } } ) # end of sapply
                       } else {
        if (p==0) { probability[1] <- 1 }
        if (p==1) { probability[m] <- 1 } } # end if ((p>0) & (q>0)) 

   output <- list(probability=probability)
   return(output) }
