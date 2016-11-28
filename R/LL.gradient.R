LL.gradient <-
function(parameter,model.type,model.name,link,ntrials,nsuccess,
                     covariates.matrix.p,covariates.matrix.scalef,
                     offset.p,offset.scalef,weights,grad.method) {
   if (grad.method=="simple") { method.args=list(eps=0.001)
                       } else { method.args=list(r=6,d=0.001) }
   gradient <- grad(LL.Regression.Binary,x=parameter,
                model.type=model.type,model.name=model.name,
                link=link,ntrials=ntrials,nsuccess=nsuccess,
                covariates.matrix.p=covariates.matrix.p,
                covariates.matrix.scalef=covariates.matrix.scalef,
                offset.p=offset.p,offset.scalef=offset.scalef,weights=weights,
                method=grad.method,method.args=method.args)
   return(gradient) }
