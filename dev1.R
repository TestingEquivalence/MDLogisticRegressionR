tmp  <- function(coef){
  a <- coef[1]
  b <- coef[2]
  a +b*tdf$x
}

distance<-function(coef){
  md$coefficients=coef
  l=predict.lm(md,df)
  logistic(l)
}
