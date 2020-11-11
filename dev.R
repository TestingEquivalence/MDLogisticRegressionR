source("dataSets.R")
source("dev1.R")
source("mdLogitRegression.R")
source("asymptoticTest.R")

#calculate p, which is  the proportion of using
df=FujiFerilitySurvey()
df$n=df$notUsing+df$using
df$p=df$using/df$n

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

#or second variant
lr1 <- glm(cbind(using,notUsing) ~ wantsMore+education+age,
          data = df, family = binomial("logit"))

mdr=nls(p ~ distance(coef), data=df,start=list(coef=lr$coefficients))

v=lr$fitted.values-df$p
sum(v*v)

w=residuals(mdr)
sum(w*w)

ww=fitted(mdr)-df$p
sqrt(sum(ww*ww))

mdl=min_dst_logit(df, "notUsing","using","wantsMore+education+age",alpha = 0.05,
                  test=asymptotic)

#checl min distance
mdl$min.distance

coef(mdr)
mdl$coef
# estimation works fine
coef(lr)

fitted(mdl)
mdl$dependentVariable

mdl$min.distance
mdl$min.epsilon
