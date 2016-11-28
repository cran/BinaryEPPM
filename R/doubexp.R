doubexp <-
function()
{ 
   linkfun <- function(mu) { leng <- length(mu)
                 sapply(1:leng, function(i) {
                    if (mu[i]<0.5) { log(2*mu[i]) 
                            } else { abs(log(2*(1-mu[i]))) }} ) }
   linkinv <- function(eta) { (1+sign(eta))/2 - sign(eta)*exp(-abs(eta))/2 }
   mu.eta <- function(eta) { exp(-abs(eta))/2 }
   valideta <- function(eta) TRUE
   link <- paste0("doubexp(",")")
   structure (list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
