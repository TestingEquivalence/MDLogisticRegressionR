

asymptStDev<-function(mdr){
  N=sum(mdr$n)
  w=mdr$n/N
  p=mdr$dependentVariable
  v=p*(1-p)/w
  dd=2*(p-fitted(mdr))
  vol=sum(dd*v*dd)/N
  return(sqrt(vol))
}

asymptoticTest<-function(mdr){
  #calculate asymptotic min eps
  vol = asymptStDev(mdr)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
