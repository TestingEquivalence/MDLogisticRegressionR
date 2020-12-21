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
  md= lm(frm, df)
  y=all.vars(as.formula(frm))[1]
  
  # logistic model for given parameters
  distance<-function(coef){
    md$coefficients=coef
    l=predict.lm(md,data)
    logistic(l)-df[[y]]
  }
  
  # calculate minimum distance estimator
  
  res=nls.lm(par=lr$coefficients, fn=distance)
  mdr$result.nls.lm=res
  
  # calculate min distance
  mdr$min.distance=sqrt(deviance(res))
  mdr$coefficients=coef(res)
  mdr$residuals=res$fvec
  mdr$fitted=res$fvec+df[[y]]
  
  # test results
  mdr$min.epsilon=NA
  
  if (asymptotic==test) {
    mdr$min.epsilon=asymptoticTest(mdr=mdr)
  }

  if ("bootstrap"==test){
    mdr$min.epsilon=bootstrapTest(df,mdr,nSimulation)
  }
  
  return(mdr)
}

updateMinDistanceModel<-function(p,mdr,df){
  n=df[[mdr$zeroCounts]]+df[[mdr$oneCounts]]
  df[[mdr$oneCounts]]=n*p
  df[[mdr$zeroCounts]]=n-df[[mdr$oneCounts]]
  
  nlr=min_dst_logit(df,zeroCounts =mdr$zeroCounts ,oneCounts =mdr$oneCounts ,
                    covariates =mdr$covariates ,test =mdr$test ,alpha =mdr$alpha )
  return(nlr)
}
