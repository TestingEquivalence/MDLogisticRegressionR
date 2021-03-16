library(forcats)

# development esoph
#########################################
df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)

df_age=aggregate(cbind(ncases,ncontrols) ~ agegp, df, sum)
write.csv(df_age,"df_age.csv")

df_alc=aggregate(cbind(ncases,ncontrols) ~ alcgp, df, sum)
write.csv(df_alc,"alcgp.csv")

df_tob=aggregate(cbind(ncases,ncontrols) ~ tobgp, df, sum)
write.csv(df_tob, "tobgp.csv")

df$alcgp=fct_collapse(df$alcgp,over80=c("80-119","120+"))
df$agegp=fct_collapse(df$agegp,i25_44=c("25-34","35-44"), over65=c("65-74","75+"))
df$tobgp=fct_collapse(df$tobgp,under10=c("0-9g/day"), over10=c("20-29","30+","10-19"))


df=aggregate(cbind(ncases,ncontrols) ~ alcgp+tobgp+agegp, df, sum)
df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ agegp+alcgp+tobgp"
lr <- glm(frm, df, family = binomial("logit"), weights = n)
lr
