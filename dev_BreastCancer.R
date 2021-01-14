################################################################
library(mlbench)
data("BreastCancer")

df=BreastCancer
df$Id=NULL
summary(df)
frm="Class~."
lm=glm(frm, df, family = binomial("logit"))
