bootstrapVolatility<-function(df,mdr,nSimulation){
  #calculate bootstrap volatility
  n=df[[mdr$zeroCounts]]+df[[mdr$oneCounts]]
  p=df[[mdr$oneCounts]]/n
  
  mdr$test="asymptotic"
  f<-function(i){
    np=resample.p(i,n,p)
    nmdr=updateMinDistanceModel(p=np,lr=mdr,df=df)
    return(nmdr$min.distance^2)
  }
  
  res=sapply(c(1:nSimulation),f)
  return(sd(res))
}

bootstrapTest1<-function(df,mdr,nSimulation){
  
  
  #calculate asymptotic min eps
  vol = bootstrapVolatility(df,mdr,nSimulation)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
bootstrapTest2proto<-function(df,mdr,nSimulation,eps){
  #find boundary point
  mdr$test="asymptotic"
  n=df[[mdr$zeroCounts]]+df[[mdr$oneCounts]]
  p=df[[mdr$oneCounts]]/n
  bp=linearBoundaryPoint(p,df,mdr,eps)
  
  #nmdr=updateMinDistanceModel(p=bp,lr=mdr,df=df)
  
  #simulate bootstrap distribution
  f<-function(i){
    np=resample.p(i,n,bp)
    nmdr=updateMinDistanceModel(p=np,lr=mdr,df=df)
    return(nmdr$min.distance)
  }
  
  res=sapply(c(1:nSimulation),f)
  
  #compute p-value
  pvalue=sum(mdr$min.distance>=res)/nSimulation
  return(pvalue)
}
bootstrapTest2<-function(df,mdr,nSimulation,lower,upper){
  f<-function(eps){
    set.seed(10071977)
    pval=bootstrapTest2proto(df,mdr,nSimulation,eps)
    return(pval-mdr$alpha)
  }
  
  res=uniroot(f,c(lower,upper),tol=mdr$alpha/100)
  return(res$root) 
}