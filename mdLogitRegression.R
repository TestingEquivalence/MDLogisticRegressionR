source("asymptoticTest.R")

# counts2logits<-function(df,zeroCounts,oneCounts,covariates){
#   fr=paste("cbind(",zeroCounts,",",oneCounts,") ~ ",covariates)
#   fr=as.formula(fr)
#   res=xtabs(fr, data = df)
#   res=ftable(res)
#   res=data.frame(expand.grid(rev(attr(res, "row.vars"))), unclass(res))
#   names(res)[names(res)=="X1"]=zeroCounts
#   names(res)[names(res)=="X2"]=oneCounts
#   return(res)
# }

logit=qlogis
logistic = plogis
asymptotic="asymptotic"
bootstrap="bootstrap"

min_dst_logit<-function(df,zeroCounts,oneCounts,covariates,  alpha=0.05,
                        test, nSimulation=1000, eps=NA){
  # prepare data
  # look for name of probability
  i=1
  p="p"
  repeat{
    if(!(p %in% colnames(df))) {break()}
    p=paste("p",i,sep="")
    i=i+1
  }
   
  kdf=data.frame(n=df[[zeroCounts]]+df[[oneCounts]])
  kdf$p=df[[oneCounts]]/kdf$n
  df[[p]]=kdf$p
  
  # formula for regression
  fr=paste(p," ~ ",covariates)
  fr=as.formula(fr)
  
  # dummy model for technical reasons
  md= lm(fr, data=df)
  
  # logistic model for given parameters
  distance<-function(coef){
    md$coefficients=coef
    l=predict.lm(md,df)
    logistic(l)
  }
  
  #logit regression for initial values
  lr = glm(fr, data = df, family = binomial("logit"), weights = kdf$n)
  
  # calculate minimum distance estimator
  fr1=paste(p,"~ distance(coef)")
  fr1=as.formula(fr1)
  mdr=nls(fr1, data=df,start=list(coef=lr$coefficients))
  
  # additional information
  mdr$zeroCounts=zeroCounts
  mdr$oneCounts=oneCounts
  mdr$test=test
  mdr$alpha=alpha
  mdr$formula=fr1
  mdr$covariates=covariates
  mdr$dependentVariable=kdf$p
  mdr$eps=eps
  mdr$n=kdf$n
  
  # calculate min distance
  mdr$min.distance=sqrt(deviance(mdr))
  mdr$coef=coef(mdr)
  names(mdr$coef)=names(coef(lr))
  
  # test results
  mdr$alpha=alpha
  mdr$min.epsilon=NA
  
  if (asymptotic==test) {
    mdr$min.epsilon=asymptoticTest(mdr=mdr)
  }

  # if ("bootstrap1"==test){
  #   mdr$min.epsilon=bootstrapTest(df,mdr,nSimulation)
  # }
  
  return(mdr)
}


