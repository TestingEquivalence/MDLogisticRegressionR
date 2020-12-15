source("asymptoticTest.R")
source("simulation.R")

simulatePowerAtModel<-function(df,n,p,lr, updateLR, nSimulation){
  set.seed(01032020)
  psim=lapply(c(1:nSimulation), resample.p,n,p)
  res=list()
  
  i=1
  repeat{
    skip= FALSE
    
    tryCatch({
      nlr=updateLR(p=psim[[i]],lr,df)
    }, 
    error = function(e) { 
      skip=TRUE
    })
    
    if(skip) { 
      next 
    }  
    
    
    res[[i]]=mdr2results(nlr)
    i=i+1
    if (i>nSimulation) break
  }
  
  return(res)
}
