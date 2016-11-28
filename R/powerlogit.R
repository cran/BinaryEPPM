powerlogit <-
function(power = 1)
{ 
   linkfun <- function(mu) { wkv <- exp(log(mu)/power)
                             log(wkv) - log(1-wkv) }
   linkinv <- function(eta) { 1 / ((1+exp(-eta)))^power }
   mu.eta <- function(eta) { power*exp(-eta) / ((1+exp(-eta)))^(power+1) }
   valideta <- function(eta) TRUE
   link <- paste0("powerlogit(", power, ")")
   structure (list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
