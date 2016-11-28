cooks.distance.BinaryEPPM <-
function (model, ...)  {
    h <- hatvalues.BinaryEPPM(model)
    k <- length(model$coefficients$p.est)
    res <- residuals.BinaryEPPM(model, type = "pearson")
    h * (res^2)/(k * (1 - h)^2) }
