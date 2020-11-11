library(tidyverse)
require(forcats)

source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")
source("bootstrapTest.R")

# data 

df=as.data.frame(Titanic)
df$nSurvived=ifelse(df$Survived=="Yes", df$Freq, 0)
df$nDeceased=ifelse(df$Survived=="No", df$Freq, 0)
df=aggregate(cbind(nSurvived,nDeceased) ~ Class+Sex+Age, df, sum)
df$n=df$nSurvived+df$nDeceased
df=df[df$n>0,]
df$p=df$nSurvived/df$n

# delete the classes child in 1st and 2nd classes, 
# because all those children survived and hence logit can not be calculated
df=df[df$p<1,]

# fitting the model and perform a single equivalence tests
###########################################################

# using logit refression
lr <- glm(p ~ Class+Sex+Age,
          data = df, family = binomial("logit"), weights = n)
lr$min.distance=n2(logit(df$p)-logit(lr$fitted.values))
write.result(lr,"lr.csv")

# using minimum distance regression
set.seed(09052020)
mdr = min_dst_logit(df,"nDeceased","nSurvived"," Class+Sex+Age",
                    test = "bootstrap1", nSimulation = 1000)
write.result(mdr,"mdr.csv")

# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities
lr = glm(p ~ Class+Sex+Age,
         data = df, family = binomial("logit"), weights = n)

mdr = min_dst_logit(df,"nDeceased","nSurvived"," Class+Sex+Age",
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


