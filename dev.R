set.seed(01012021)
rp=randomExteriorPoint(p=df$proportionUsing,df,mdr,updateLR = updateMinDistanceModel,eps = 0.4)
nmdr=updateMinDistanceModel(rp,mdr,df)
nmdr$min.distance
