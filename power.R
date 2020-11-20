randomExteriorPoint<-function(p,df,mdr, eps){
  repeat{
    skip= FALSE
    sp=resample.p(1,mdr$n,p)
    
    tryCatch({mdr= updateMinDistanceModel(sp,mdr,df)}, 
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


simulatePowerAtBoundary<-function(df,p,mdr, nSimulation, eps){
  set.seed(01032020)
  exteriorPoints=list()
  bndPoints=list()
  test=mdr$test
  mdr$test=asymptotic
  nPoints=100
  
  for (i in c(1:nPoints)){
    exteriorPoints[[i]]=randomExteriorPoint(p,df,mdr,eps)
  }
  
  for (i in c(1:nPoints)){
    bndPoints[[i]]=linearBoundaryPoint(interiorPoint = p,
                                       exteriorPoint = exteriorPoints[[i]],
                                       df, mdr, eps)
  }
  
  mdr$test=test
  
  cl=getCluster()
  power=parSapply(cl,bndPoints, simulatePowerAtPoint,df,mdr, nSimulation,eps)
  stopCluster(cl)
  
  # power=sapply(bndPoints, simulatePowerAtPoint,df,mdr, nSimulation,eps)
  
  return(power)
  }
  
  

