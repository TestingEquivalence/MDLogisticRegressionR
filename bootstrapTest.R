bootstrapVolatility<-function(df,mdr,nSimulation){
  #calculate bootstrap volatility
  n=df[[mdr$zeroCounts]]+df[[mdr$oneCounts]]
  p=df[[mdr$oneCounts]]/n
  
  mdr$test="asymptotic"
  res=rep(NA,nSimulation)
  i=1
  
  repeat{
    skip= FALSE
    np=resample.p(i,n,p)
    
    tryCatch({
      nmdr=updateMinDistanceModel(p=np,mdr=mdr,df=df)
      }, 
    error = function(e) { 
      skip=TRUE
      })
    
    if(skip) { 
      next 
    }  
    
    res[i]=nmdr$min.distance^2
    i=i+1
    if (i>nSimulation) break
  }
  
  
  
  return(sd(res))
}

bootstrapTest<-function(df,mdr,nSimulation){
  
  
  #calculate asymptotic min eps
  vol = bootstrapVolatility(df,mdr,nSimulation)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
