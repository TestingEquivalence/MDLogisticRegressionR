require(tidyverse)
require(forcats)


source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")
source("bootstrapTest.R")

# prepare data
df=esoph
#View(df)
df=aggregate(cbind(ncases,ncontrols) ~ agegp+alcgp, df, sum)
df$agegp=factor(df$agegp, ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
str(df)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

formula_lr="p ~ agegp+alcgp"
formula_mdr="agegp+alcgp"

# fitting the model and perform a single equivalence tests
###########################################################

# using logit refression
lr <- glm(formula_lr,
          data = df, family = binomial("logit"), weights = n)
lr$min.distance=sqrt(sum((df$p-lr$fitted.values)^2))
write.result(lr,"lr.csv")

# using minimum distance regression
set.seed(01012021)
mdr = min_dst_logit(df,"ncontrols","ncases",formula_mdr,
                    test = asymptotic, nSimulation = 200)
write.result(mdr,"mdr.csv")

# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities
lr = glm(formula_lr,
         data = df, family = binomial("logit"), weights = n)

mdr = min_dst_logit(df,"ncontrols","ncases",formula_mdr,
                    test = asymptotic, nSimulation = 200)

# compute distribution using logit regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=fitted(lr),
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_lr_power_lr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=fitted(mdr),
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_lr.csv")

# compute distribution using minimum distance regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=fitted(lr),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_lr_power_mdr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=fitted(mdr),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_mdr.csv")

# compute test power at the fitted model 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(df,"ncontrols","ncases",formula_mdr,
                    test = asymptotic, nSimulation = 200)

res=simulatePowerAtModel(df,
                         n=df$n,
                         p=fitted(mdr),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"size_mdr.csv")

# compute test power at the random boundary points 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(df,"ncontrols","ncases",formula_mdr,
                    test = asymptotic, nSimulation = 200)

res= simulatePowerAtBoundary(df,p=df$p,mdr, nSimulation=1000, eps=0.6)
write.csv(res,"power_mdr.csv")

