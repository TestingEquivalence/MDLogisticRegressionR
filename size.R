source("asymptoticTest.R")
source("simulation.R")

simulatePowerAtModel<-function(df,n,p,lr, updateLR, nSimulation){
  set.seed(01032020)
  psim=lapply(c(1:nSimulation), resample.p,n,p)
  
  f<-function(p,lr){
    nlr=updateLR(p,lr,df)
    res=mdr2results(nlr)
    return(res)
  }
  
  lres=lapply(psim,f,lr)
  return(lres)
}
