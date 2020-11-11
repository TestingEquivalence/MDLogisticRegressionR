randomExteriorPoint<-function(df,n,p,lr, updateLR, eps){
  repeat{
    sp=resample.p(1,n,p)
    nlr=updateLR(sp,lr,df)
    if (nlr$min.distance>=eps){
      return(sp)
    }
  }
}

linearBoundaryPoint<-function(exteriorPoint,df,mdr, eps){
   nmdr=updateMinDistanceModel(p=exteriorPoint,lr=mdr,df=df)
  a=eps/nmdr$min.distance
  la=a*logit(exteriorPoint)+(1-a)*nmdr$fitted.values
  pa=logistic(la)
  #pmdr=updateMinDistanceModel(p=pa,lr=mdr,df=df)
  return(pa)
}

simulatePowerAtBoundary<-function(df,n,p,mdr, nSimulation, eps){
  set.seed(01032020)
  exteriorPoints=list()
  test=mdr$test
  mdr$test=""
  
  for (i in c(1:100)){
    exteriorPoints[[i]]=randomExteriorPoint(df,n,p,lr=mdr, updateLR=updateMinDistanceModel, eps)
  }
  
  mdr$test=test
  bndPoints=lapply(exteriorPoints, linearBoundaryPoint, df, mdr, eps)
  
  cl=getCluster()
  power=parSapply(cl,bndPoints, simulatePowerAtPoint,df,n,mdr, nSimulation,eps)
  stopCluster(cl)
  
  # power=sapply(bndPoints, simulatePowerAtPoint,df,n,mdr, nSimulation,eps)
  
  return(power)
  }
  
  

