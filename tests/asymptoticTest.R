

asymptStDev<-function(mdr){
  N=sum(mdr$weights)
  w=mdr$weights/N
  y=all.vars(as.formula(mdr$frm))[1]
  p=mdr$data[[y]]
  v=p*(1-p)/w
  dd=2*mdr$residuals
  vol=sum(dd*v*dd)
  return(sqrt(vol))
}

asymptoticTest<-function(mdr){
  #calculate asymptotic min eps
  N=sum(mdr$weights)
  vol = asymptStDev(mdr)/sqrt(N)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
