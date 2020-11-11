bootstrapVolatility<-function(df,mdr,nSimulation){
  #calculate bootstrap volatility
  n=df[[mdr$zeroCounts]]+df[[mdr$oneCounts]]
  p=df[[mdr$oneCounts]]/n
  
  mdr$test="asymptotic"
  f<-function(i){
    np=resample.p(i,n,p)
    nmdr=updateMinDistanceModel(p=np,mdr=mdr,df=df)
    return(nmdr$min.distance^2)
  }
  
  res=sapply(c(1:nSimulation),f)
  return(sd(res))
}

bootstrapTest<-function(df,mdr,nSimulation){
  
  
  #calculate asymptotic min eps
  vol = bootstrapVolatility(df,mdr,nSimulation)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
