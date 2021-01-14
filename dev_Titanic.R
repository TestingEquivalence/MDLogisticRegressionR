# development titanic
#########################################

df=Titanic
frm="Survived ~ Class+Sex+Age"
lr <- glm(frm, df, family = binomial("logit"), weights = Freq)

summary(lr)
confint(lr)

df=as.data.frame(Titanic)
df$Survived=ifelse(df$Survived=="Yes", df$Freq, 0)
df$Deceased=ifelse(df$Survived==0, df$Freq, 0)
df=aggregate(cbind(Survived,Deceased) ~ Class+Sex+Age, df, sum)
df$n=df$Survived+df$Deceased
df$p=df$Survived/df$n

frm="p ~ Class+Sex+Age"
lr1 <- glm(frm, df, family = binomial("logit"), weights = n)
summary(lr1)
confint(lr1)

df10=df[df$n>=10,]
lr10 <- glm(frm, df10, family = binomial("logit"), weights = n)
summary(lr10)
confint(lr10)

df20=df[df$n>=20,]
lr20 <- glm(frm, df20, family = binomial("logit"), weights = n)
summary(lr20)
confint(lr20)

dfAdO=df[df$Age=="Adult",]
dfAdO$Age=NULL
frm="p ~ Class+Sex"
lrAdO <- glm(frm, dfAdO, family = binomial("logit"), weights = n)
summary(lrAdO)
confint(lrAdO)


