logit<-qlogis
logistic <- plogis

n2<-function(v){
  v=v*v
  return(sqrt(sum(v)))
}

l2<-function(v1,v2){
  return(n2(v1-v2))
}

ll2<-function(p,q){
  l2(logit(p),logit(q))
}

asymptStDev<-function(df,mdr){
  N=sum(df$n)
  w=df$n/N
  v=w*df$p*(1-df$p)
  v=1/v
  v=v*4*(df$l-mdr$fitted.values)*(df$l-mdr$fitted.values)
  vol=sum(v)/N
  return(sqrt(vol))
}

asymptoticTest<-function(df,mdr){
  #calculate asymptotic min eps
  vol = asymptStDev(df,mdr)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
