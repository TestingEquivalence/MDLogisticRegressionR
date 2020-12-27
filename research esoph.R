source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")
source("bootstrapTest.R")

# prepare data
df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ tobgp+alcgp, df, sum)
str(df)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)
#df=df[df$n>=10,]

frm="p ~ tobgp+alcgp"

# fitting the model and perform a single equivalence tests
###########################################################

# using logit refression
lr <- glm(frm, df, family = binomial("logit"), weights = n)
lr$min.distance=sqrt(sum((df$p-lr$fitted.values)^2))
write.result(lr,"lr.csv")

# using minimum distance regression
set.seed(01012021)
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)
write.result(mdr,"mdr.csv")



# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities

lr = glm(frm,df, family = binomial("logit"), weights =n)

mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

# compute distribution using logit regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=lr$fitted.values,
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_lr_power_lr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=mdr$fitted,
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
                         p=mdr$fitted,
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_mdr.csv")

# compute test power at the fitted model 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

res=simulatePowerAtModel(df,
                         n=df$n,
                         p=fitted(mdr),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"size_mdr.csv")

# compute test power at the random boundary points 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

res= simulatePowerAtBoundary(p=df$p,mdr, nSimulation=1000, eps=0.34)
write.csv(res,"power_mdr.csv")
