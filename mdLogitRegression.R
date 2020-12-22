source("asymptoticTest.R")
source("bootstrapTest.R")
source("simulation.R")
library(minpack.lm)

logit=qlogis
logistic = plogis
asymptotic="asymptotic"
bootstrap="bootstrap"

min_dst_logit<-function(formula,data, weights,  test, alpha=0.05,
                        nSimulation=200){
  
  #initial information
  mdr=list()
  mdr$data=data
  mdr$formula=as.formula(formula)
  mdr$frm=formula
  mdr$weights=weights
  mdr$alpha=alpha
  mdr$test=test
  mdr$nSimulation=nSimulation
  
 
  
  
  
  #logit regression for initial values
  lr <- glm(mdr$formula,mdr$data, family = quasibinomial("logit"), weights =mdr$weights)
  
  # dummy model for technical reasons
  md= lm(mdr$frm, mdr$data)
  y=all.vars(as.formula(mdr$frm))[1]
  
  # logistic model for given parameters
  distance<-function(coef){
    md$coefficients=coef
    l=predict.lm(md,mdr$data)
    logistic(l)-mdr$data[[y]]
  }
  
  # calculate minimum distance estimator
  
  res=nls.lm(par=lr$coefficients, fn=distance)
  mdr$result.nls.lm=res
  
  # calculate min distance
  mdr$min.distance=sqrt(deviance(res))
  mdr$coefficients=coef(res)
  mdr$residuals=res$fvec
  mdr$fitted=res$fvec+mdr$data[[y]]
  
  # test results
  mdr$min.epsilon=NA
  
  if (asymptotic==test) {
    mdr$min.epsilon=asymptoticTest(mdr=mdr)
  }

  if ("bootstrap"==test){
    mdr$min.epsilon=bootstrapTest(mdr,nSimulation)
  }
  
  return(mdr)
}

updateMinDistanceModel<-function(p,mdr){
  df=mdr$data
  y=all.vars(as.formula(mdr$frm))[1]
  df[[y]]=p
  
  nlr=min_dst_logit(mdr$frm,df,weights=mdr$weights,test = mdr$test)
  return(nlr)
}
