source("dataSets.R")
source("mdLogitRegression.R")

# This project provides the minimum distance logistic regression for the situation
# that the response is binary (independent variable is binary)
# and all covariates (all dependent variables) are categorical.

# The minimum distance logistic regression uses the Euclidean distance between
# the observed counting frequencies p(n,x) and the corresponding logistic regression probabilities p(b,x),
# where n is the number of observation, x is a value of covariates and
# b is the vector of coefficients.
# The regression coefficients b are estimated using the minimization of this distance.

# In addition, this project provides the means to demonstrate that 
# the fitted logistic model is sufficiently close to the true underlying distribution.

# Let us consider the equivalence test problem:
# H0={ distance between the model und true underlying distrubition is larger than epsilon}
# H1={ distance between the model und true underlying distrubition is smaller than epsilon},
# where epsilon is the tolerance parameter.
# If H0 can be rejected for an appropriate value of epsilon, then the logistic regression model
# is sufficiently close to the true underlying distribution. 
# This program computes the minimum tolerance parameter epsilon, for which H0 can be rejected
# and hence the equivalence can be established.
# There are 4 different methods to compute the minimum tolerance parameter epsilon at the given 
# significance level:
# - approximation by asymptotic distribution, uses also asymptotic variance
# - approximation by asymptotic distribution using bootstrapped variance
# - empirical bootstrap method
# - studentized bootstrap method, also known as bootstrap-t method



# -----------------------------------------------------------
# Example Fiji Fertility Study
# -----------------------------------------------------------

# See dataSets.R for the detailed information on this data set
# Prepare data:

df=FujiFerilitySurvey()
str(df)

df$n=df$using+df$notUsing
df$p=df$using/df$n
frm="p ~ wantsMore+education+age"

# Perform the usual logistic regression:
lr = glm(frm,df, family = binomial("logit"), weights =n)

# Compute the distance between the observed counting frequencies 
# and the model, computed by usual logistic regression: 
lr$min.distance=sqrt(sum((df$p-lr$fitted.values)^2))

# Perform the minimum distance logistic regression:
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

# Function min_dst_logit has the following arguments:
# formula - an object of class "formula": a symbolic description of the model to be fitted
# data - data frame containing the variables in the model
# weights - the vector of all counts for each combinations of predictors, i.e. number of observations in each cell of multi-way contingency table 
# test - the type of the equivalence test to be performed, possible values are: 
#     asymptotic; asymptoticBootstrapVariance; empiricalBootstrap; tPercentileBootstrap
# alpha - the nominal significance level for the equivalence test
# nSimulation - number of the bootstrap simulations for the equivalence test

# Compare the regression coefficients, which are estimated by the usual logistic regression 
# and by the minimum distance logistic regression:
coef(lr)
coef(mdr)

# Compare the Euclidean distances between the models and data:
lr$min.distance
mdr$min.distance

# Perform the equivalence test using all available methods and compare the results
mdr_asympt= min_dst_logit(frm,df,weights=df$n,test = asymptotic)
mdr_asympt_bootstrap_variance=min_dst_logit(frm,df,weights=df$n,test = asymptoticBootstrapVariance, nSimulation = 1000)
mdr_empirical_bootstrap=min_dst_logit(frm,df,weights=df$n,test = empiricalBootstrap, nSimulation = 1000)
mdr_tPercentile_bootstrap=min_dst_logit(frm,df,weights=df$n,test = tPercentileBootstrap, nSimulation = 1000)

# Compare the results. The minimum tolerance parameter epsilon, for which H0 can be rejected, is reported.
mdr_asympt$min.epsilon
mdr_asympt_bootstrap_variance$min.epsilon
mdr_empirical_bootstrap$min.epsilon
mdr_tPercentile_bootstrap$min.epsilon
