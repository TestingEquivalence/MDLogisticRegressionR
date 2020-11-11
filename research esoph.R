require(tidyverse)
require(forcats)

source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")
source("bootstrapTest.R")

# data
df=esoph
df$agegp=fct_collapse(df$agegp,
                       "25-44"=c("25-34","35-44"),
                       "45-64"=c("45-54","55-64"),
                       "65+"=c("65-74","75+"))
df$alcgp=fct_collapse(df$alcgp,
                       "0-79"=c("0-39g/day","40-79"),
                       "80+"=c("80-119","120+"))
df$tobgp=fct_collapse(df$tobgp,
                       "0-19"=c("0-9g/day","10-19"),
                       "20+"=c("20-29","30+"))


df=aggregate(cbind(ncases,ncontrols) ~ agegp+alcgp+tobgp, df, sum)
str(df)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

# fitting the model and perform a single equivalence tests
###########################################################

# using logit refression
lr <- glm(p ~ agegp+alcgp+tobgp,
          data = df, family = binomial("logit"), weights = n)
lr$min.distance=n2(logit(df$p)-logit(lr$fitted.values))
write.result(lr,"lr.csv")

# using minimum distance regression
set.seed(09052020)
mdr = min_dst_logit(df,"ncontrols","ncases","agegp+alcgp+tobgp",
                    test = "bootstrap2", nSimulation = 1000)
write.result(mdr,"mdr.csv")

# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities
lr = glm(p ~ agegp+alcgp+tobgp,
         data = df, family = binomial("logit"), weights = n)

mdr = min_dst_logit(df,"ncontrols","ncases","agegp+alcgp+tobgp",
                    test = "asymptotic")

# compute distribution using logit regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=lr$fitted.values,
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_lr_power_lr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=logistic(mdr$fitted.values),
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_lr.csv")

# compute distribution using minimum distance regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=lr$fitted.values,
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_lr_power_mdr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=logistic(mdr$fitted.values),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_mdr.csv")

