source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")
source("bootstrapTest.R")

# data
df=FujiFerilitySurvey()
str(df)

df$n=df$using+df$notUsing
df$proportionUsing=df$using/df$n

# fitting the model and perform a single equivalence tests
###########################################################

# using logit refression
lr <- glm(proportionUsing ~ wantsMore+education+age,
                data = df, family = binomial("logit"), weights = n)
lr$min.distance=n2(logit(df$proportionUsing)-logit(lr$fitted.values))
write.result(lr,"lr.csv")

# using minimum distance regression
set.seed(09052020)
mdr = min_dst_logit(df,"notUsing","using","wantsMore+education+age",
                         test = bootstrap2, nSimulation = 1000)
write.result(mdr,"mdr.csv")


# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities
lr = glm(proportionUsing ~ wantsMore+education+age,
         data = df, family = binomial("logit"), weights = n)

mdr = min_dst_logit(df,"notUsing","using","wantsMore+education+age",
                    test = asymptotic, nSimulation = 1000)

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

# compute test power at the fitted model 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(df,"notUsing","using","wantsMore+education+age",
                        test = bootstrap2pvalue, nSimulation = 1000, eps=2)

res=simulatePowerAtModel(df,
                         n=df$n,
                         p=logistic(mdr$fitted.values),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"size_mdr.csv")

# compute test power at the random boundary points 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(df,"notUsing","using","wantsMore+education+age",
                    test = bootstrap2pvalue, nSimulation = 1000, eps=2)

res= simulatePowerAtBoundary(df,n=df$n,
                             p=logistic(mdr$fitted.values),
                             mdr,nSimulation=1000,eps=2)
write.csv(res,"power_mdr.csv")
