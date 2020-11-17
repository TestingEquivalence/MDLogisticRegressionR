randomExteriorPoint<-function(p,df,mdr, updateLR, eps){
  repeat{
    skip= FALSE
    sp=resample.p(1,mdr$n,p)
    
    tryCatch({mdr=updateLR(sp,mdr,df)}, 
             error = function(e) { skip=TRUE})
    
    if(skip) { 
      next 
    }  
    
    
    if (mdr$min.distance>=eps){
      return(sp)
    }
  }
}

linearBoundaryPoint<-function(interiorPoint,exteriorPoint,df,mdr, eps){
  mdr$test=asymptotic
  aim<-function(a){
    #print(a)
    lc=a*interiorPoint+(1-a)*exteriorPoint
    nmdr=updateMinDistanceModel(lc,mdr,df)
    return(nmdr$min.distance-eps)
  }
  
  a=uniroot(aim, c(0,1))$root
  return(a*interiorPoint+(1-a)*exteriorPoint)
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
  
  

