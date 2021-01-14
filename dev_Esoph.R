# development esoph
#########################################
df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ alcgp+tobgp+agegp, df, sum)
str(df)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ tobgp+alcgp+agegp"
lr=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr)
confint(lr)

##################
# age only

df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ agegp, df, sum)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ agegp"
lr_age=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr_age)
confint(lr_age)


##################
# alc only

df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ alcgp, df, sum)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ alcgp"
lr_alc=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr_alc)
confint(lr_alc)

##################
# tob only

df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ tobgp, df, sum)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ tobgp"
lr_tob=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr_tob)
confint(lr_tob)

##################
# age +alc 

df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ agegp+alcgp, df, sum)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)
df=df[df$n>=20,]

frm="p ~ agegp+alcgp"
lr_age_alc=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr_age_alc)
confint(lr_age_alc)

##################
# age +tob 

df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ agegp+tobgp, df, sum)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ agegp+tobgp"
lr_age_tob=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr_age_tob)
confint(lr_age_tob)

##################
# alc +tob 

df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)
df=aggregate(cbind(ncases,ncontrols) ~ alcgp+tobgp, df, sum)

df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ alcgp+tobgp"
lr_alc_tob=glm(frm, df, family = binomial("logit"), weights = n)
summary(lr_alc_tob)
confint(lr_alc_tob)

