CBprob <-
function(twoparameter,nt) {
   p   <- twoparameter[1] 
   q   <- 1 - p 
   rho <- twoparameter[2]
   probability <- rep(0,(nt + 1))
   if ((p>0) & (q>0)) {
           probability <- sapply(0:nt, function(i) { 
# equation (4) of Kupper & Haseman, Biometrics (1978) with theta replaced by rho
# so the denominator is 2*p*q not 2*(pq)**2
# check on limits done in Model.BCBinProb
           wks <- 1 + rho*( (i-nt*p)**2 + i*(2*p-1) - nt*p*p )/(2*p*q)
           probability[i+1] <- wks*exp(i*log(p) + (nt-i)*log(q)) } ) # end of sapply
                       } else {
        if (p<=0) { probability[1] <- 1 }
        if (p>=1) { probability[nt + 1] <- 1 } } # end if ((p>0) & (q>0)) 

   output <- list(probability=probability)
   return(output) }
