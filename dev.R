source("dataSets.R")
source("dev1.R")

logit<-qlogis
logistic <- plogis

#calculate p, which is  the proportion of using
df=FujiFerilitySurvey()
df$n=df$notUsing+df$using
df$p=df$using/df$n


x <- 1:10
y <- 2*x + 3
yeps <-y+  rnorm(length(x), sd = 0.01)  # added noise
tdf=data.frame(y=yeps,x=x,z=2*x)
nls(y ~ a + b*x,data=tdf, start = list(a = 0.12345, b = 0.54321))#

nls(y ~ tmp(coef),data=tdf, start = list(coef = c(0.12345, 0.54321)))



#dummy model for technical reasons
md= lm(p ~ wantsMore+education+age, data=df)

#test of calculation tricks using manual regression coefficients
wmd=md
wmd$coefficients=rep(0,length(wmd$coefficients))
predict(md,df)
predict.lm(wmd,df)

# logit regression for initial values
lr <- glm(p ~ wantsMore+education+age,
          data = df, family = binomial("logit"), weights = n)

mdr=nls(p ~ distance(coef), data=df,start=list(coef=lr$coefficients))

v=lr$fitted.values-df$p
sum(v*v)
