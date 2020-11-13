randomExteriorPoint<-function(p,df,mdr, updateLR, eps){
  repeat{
    skip= FALSE
    sp=resample.p(1,mdr$n,p)
    
    tryCatch({mdr=updateLR(sp,mdr,df)}, 
             error = function(e) { skip=TRUE})
    
    if(skip) { 
      print("skipped")
      next 
    }  
    
    
    if (mdr$min.distance>=eps){
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

linearBoundaryPoint2<-function(interiorPoint,exteriorPoint,df,mdr, eps){
 
  aim<-function(a){
    lc=linComb(P,Q,a)
    res = nearestPowerLaw(lc,kmin,kmax,1,3)
    beta=res$minimum
    distance=res$objective
    return(distance-eps)
  }
  
  aMin=uniroot(aim, c(0,1))
  return(linComb(p,q,aMin$root))
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
  
  

