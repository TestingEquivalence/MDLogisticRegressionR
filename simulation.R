library(parallel)

resample.p<-function(i,n,p){
  f<-function(r){
    return(rbinom(1,r[1],r[2])/r[1])
  }
  
  pdf=data.frame(n=n,p=p)
  repeat{
    res=apply(X=pdf, MARGIN=1,FUN = f)
    res=as.vector(res)
    ok=TRUE
    for (r in res){
      if (r<=0) {ok=FALSE}
      if (r>=1) {ok=FALSE}
    }
    if (ok) {break}
  }
  
  return(res)
}

mdr2results<-function(mdr){
  #gather resuts for effcience purpose
  res=list(min.distance=mdr$min.distance,
           min.epsilon=mdr$min.epsilon,
           coefficients=mdr$coefficients
           )

}

getTitle<-function(res){
  #title
  title=paste("min_distance"," min_epsilon",sep=",")
  
  for (cf in names(res[[3]])) {
    title=paste(title,cf,sep=",")
  }
  
  return(title)
}

res2row<-function(row){
  sl=paste(row[[1]],row[[2]], sep=",")
  for (cf in row[[3]]){
    sl=paste(sl,cf,sep=",")
  }
  return(sl)
}

write.result<-function(mdr,fname){
  res=mdr2results(mdr)
  fc=file(fname)
  title=getTitle(res)
  row=res2row(res)
  rows=c(title,row)
  writeLines(rows,con=fc)
  close(fc)
}

write.results<-function(res,fname){
    fc=file(fname)
    title=getTitle(res[[1]])
    
    #data
    rows=sapply(res, res2row)
    rows=c(title,rows)
    #rows=as.vector(rows)
    
    writeLines(rows,con=fc)
    close(fc)
}

updateLogitModel<-function(p,lr,df){
  df=lr$data
  frm=lr$formula
  depVar=all.vars(frm)[1]
  
  df[[depVar]]=p
  nlr=update(lr,data=df)
  nlr$min.distance=n2(logit(p)-logit(nlr$fitted.values))
  if (is.infinite(nlr$min.distance)){browser()}
  
  return(nlr)
}

updateMinDistanceModel<-function(p,lr,df){
  n=df[[lr$zeroCounts]]+df[[lr$oneCounts]]
  df[[lr$oneCounts]]=n*p
  df[[lr$zeroCounts]]=n-df[[lr$oneCounts]]
  
  nlr=min_dst_logit(df,zeroCounts =lr$zeroCounts ,oneCounts =lr$oneCounts ,
                covariates =lr$covariates ,test =lr$test ,alpha =lr$alpha, eps=lr$eps )
  return(nlr)
}

simulatePowerAtPoint<-function(p,df,n,mdr, nSimulation,eps){
  set.seed(01032020)
  psim=lapply(c(1:nSimulation), resample.p,n,p)
  
  f<-function(np){
    nmdr=updateMinDistanceModel(p=np,lr=mdr,df=df)
    
    if (mdr$test=="bootstrap2pvalue"){
      return (nmdr$min.epsilon<=mdr$alpha)
    }
    
    return(nmdr$min.epsilon<=eps)
  }
  
  res=sapply(psim,f)
  return(sum(res==TRUE)/nSimulation)
}


# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("min_dst_logit","resample.p","updateMinDistanceModel","simulatePowerAtPoint",
                     "logit", "logistic","n2","l2","ll2","asymptStDev","asymptoticTest",
                     "bootstrapVolatility","bootstrapTest1","bootstrapTest2proto",
                     "bootstrapTest2","linearBoundaryPoint"))
  
  return(cl)
}