BBprob <-
function(twoparameter,nt) {
   p     <- twoparameter[1]
   q     <- 1 - p 
   theta <- twoparameter[2]
   probability <- rep(0,(nt + 1))

   if ((p>0) & (q>0)) {
# Calculating probabilities using Smith AS189 formula
# product of (1 + r*theta) r=0 to r=n-1; constant
# as a single sample of grouped binary data
      ntm1  <- nt - 1
      if (ntm1>0) { denom <- sum(log(abs(c(rep(1,ntm1))+c(1:ntm1)*c(rep(theta,ntm1)))))
          } else { denom <- 0 }
      probability <- sapply(0:nt, function(i) { 
# the abs( ) in the following expressions is because sometimes rounding errors 
# causes small numbers less than 0 for the largest value of r multiplying theta
# product of (p + r*theta) r=0 to r=x-1
         xm1 <- i - 1
         if (xm1>-1) { numer_one <- sum(log(abs(c(rep(p,i))+c(0:xm1)*c(rep(theta,i)))))
                     } else { numer_one <- 0 }
# product of (1 - p + r*theta) r=0 to r=n-x-1
         nmx   <- nt - i
         nmxm1 <- nmx - 1
         if (nmxm1>-1) { numer_two <- sum(log(abs(c(rep((1-p),nmx))+c(0:nmxm1)*c(rep(theta,nmx)))))
                       } else { numer_two <- 0 } 
      probability[i+1] <- exp(numer_one + numer_two - denom) } ) # end of sapply
                       } else {
        if (p<=0) { probability[1] <- 1 }
        if (p>=1) { probability[nt + 1] <- 1 } } # end if ((p>0) & (q>0)) 
   output <- list(probability=probability)
   return(output) }
