set.seed(01022021)
rp=randomExteriorPoint(p=df$proportionUsing,df,mdr,eps = 0.45)
nmdr=updateMinDistanceModel(rp,mdr,df)
nmdr$min.distance

lrp=linearBoundaryPoint(df$proportionUsing,rp,df,mdr, 0.45)
nmdr=updateMinDistanceModel(lrp,mdr,df)
nmdr$min.distance
