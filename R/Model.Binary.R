Model.Binary <-
function(parameter,model.type,model.name,link,ntrials,
                   covariates.matrix.p,covariates.matrix.scalef,
                   offset.p,offset.scalef) {
   if (model.type=="p only") { 
      if ((model.name=="binomial") | (model.name=="generalized binomial")) { 
            output <- Model.GB(parameter,model.name,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               offset.p=offset.p) 
                                         } # end of binomial, generalized binomial
      if ((model.name=="beta binomial") | (model.name=="correlated binomial")) { 
            output <- Model.BCBinProb(parameter,model.type,model.name,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               offset.p=offset.p)
                                         } # end of model BB or CB
                                } # end of p only models
   if (model.type=="p and scale-factor") { 
      if (model.name=="generalized binomial") { 
            output <- Model.JMVGB(parameter,model.name,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               covariates.matrix.scalef=covariates.matrix.scalef,
                               offset.p=offset.p,offset.scalef=offset.scalef) 
                                         } # end of model.JMVGB
      if ((model.name=="beta binomial") | (model.name=="correlated binomial")) { 
            output <- Model.BCBinProb(parameter,model.type,model.name,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               covariates.matrix.scalef=covariates.matrix.scalef,
                               offset.p=offset.p,offset.scalef=offset.scalef) 
                                         } # end of model BB or CB
                                } # end of p and scale-factor models
   return(output) }
