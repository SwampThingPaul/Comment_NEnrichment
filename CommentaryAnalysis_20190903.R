## 
## Code to support data/information specific 
## to commentary {title}
##
##
## Code was compiled by Paul Julian on 2019-09-03
## contact info: pjulian@ufl.edu

#Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(wesanderson)

# MDL Handling ------------------------------------------------------------
N.mw=14.0067

dates<-as.Date(c("1978-05-01","2019-05-01"))
dat<-DBHYDRO_WQ(dates[1],dates[2],"P36",18)
dat<-subset(dat,Collection.Method=="G");#Grab sample only
dat$WY<-WY(dat$Date.EST);#Water Year
dat$MDL.rm<-with(dat,ifelse(Value<0,NA,Value));#values <MDL replaced with NA

ann.mean<-ddply(dat,"WY",summarise,mean.val.HalfMDL=mean(HalfMDL*N.mw,na.rm=T),N.val.HalfMDL=N(HalfMDL),mean.val.rm=mean(MDL.rm*N.mw,na.rm=T),N.val.rm=N(MDL.rm))

#Trend analysis
with(ann.mean,cor.test(mean.val.rm,WY,method="kendall",exact=F))
with(ann.mean,cor.test(mean.val.HalfMDL,WY,method="kendall",exact=F))

#Plot to visualize information
#Figure 1
par(family="serif")
plot(mean.val.HalfMDL~WY,ann.mean,ylim=c(0,1.25),xlim=c(1985,2019),type="n",ylab="Nitrate-Nitrite (mM)",xlab="Water Year",las=1)
with(ann.mean,lines(WY,mean.val.rm,lty=2,col="grey",lwd=2))
with(ann.mean,points(WY,mean.val.rm,pch=21,bg="grey",cex=1.25))
with(ann.mean,lines(WY,mean.val.HalfMDL,lty=2,col=adjustcolor("dodgerblue1",0.25),lwd=2))
with(ann.mean,points(WY,mean.val.HalfMDL,pch=21,bg=adjustcolor("dodgerblue1",0.5),cex=1.25))

leg.text=c("< MDL replaced with \u00BD MDL","< MDL data removed")
legend("topright",legend=leg.text,pch=c(NA,NA),col=c("grey",adjustcolor("dodgerblue1",0.5)),pt.bg=c(NA,NA),lwd=c(1,1),lty=c(2,2),pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,text.col="white")
legend("topright",legend=leg.text,pch=c(21,21),col=c("black"),pt.bg=c("grey",adjustcolor("dodgerblue1",0.25)),lwd=c(0.1,0.1),lty=c(NA,NA),pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)


# Boxplot-outlier ---------------------------------------------------------
set.seed(1)
#Synthetic dataset
testIQR<-data.frame()
for(i in seq(1,5000,2)){
  tmp<-rnorm(i)
  q.val<-quantile(tmp,c(0.25,0.5,0.75))
  IQR<-as.numeric((q.val[3]-q.val[1]))
  N.out<-N(boxplot(tmp)$out)
  testIQR<-rbind(testIQR,data.frame(samp=i,IQR.val=IQR,N.out=N.out))
}


par(family="serif",mar=c(1,3.5,0.5,0.5),oma=c(2.5,1,0.25,0.75));
layout(matrix(c(1:2),2,1,byrow=T))

xlim.val<-c(0,5000);by.x<-2000;xmaj<-seq(xlim.val[1],xlim.val[2],by.x);xmin<-seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val<-c(0,2);by.y<-0.5;ymaj<-seq(ylim.val[1],ylim.val[2],by.y);ymin<-seq(ylim.val[1],ylim.val[2],by.y/2)
plot(IQR.val~samp,testIQR,type="n",ylim=ylim.val,xlim=xlim.val,ylab=NA,xlab=NA,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(testIQR,points(samp,IQR.val,pch=21,col=adjustcolor("indianred1",0.5),bg=adjustcolor("indianred1",0.25),lwd=0.5,cex=0.75))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Inter-quartile Range")

xlim.val<-c(0,5000);by.x<-2000;xmaj<-seq(xlim.val[1],xlim.val[2],by.x);xmin<-seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val<-c(0,60);by.y<-20;ymaj<-seq(ylim.val[1],ylim.val[2],by.y);ymin<-seq(ylim.val[1],ylim.val[2],by.y/2)
plot(N.out~samp,testIQR,type="n",ylim=ylim.val,xlim=xlim.val,ylab=NA,xlab=NA,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(testIQR,points(samp,N.out,pch=21,col=adjustcolor("dodgerblue1",0.5),bg=adjustcolor("dodgerblue1",0.25),lwd=0.5,cex=0.75))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.75,"Number of Outliers")
mtext(side=1,line=1.75,"Sample Size")


# Rainfall ----------------------------------------------------------------
rf.sites<-data.frame(SITE=c("BEELNE","BEELNE","S5A","L005","3A-S_R","3A-S_R","NP-P36","NP-P37"),
                    DBKEY=c("05963","TY244","15202","12515","05865","HC941","06038","H2001"),
                    Region=c("Orlando","Orlando","West Palm","Lake Okeechobee","Southern WCA3A","Southern WCA3A","ENP (SRS)", "ENP (Taylor Slough)"))

dates<-as.Date(c("1983-05-01","2019-05-01"))
rf.dat<-data.frame()
for(i in 1:nrow(rf.sites)){
  tmp<-DBHYDRO_daily(dates[1],dates[2],rf.sites$DBKEY[i])
  tmp$DBKEY<-as.character(rf.sites$DBKEY[i])
  rf.dat<-rbind(tmp,rf.dat)
  print(i)
}
rf.dat<-merge(rf.dat,rf.sites,"DBKEY")
rf.dat2<-ddply(rf.dat,c("SITE","Region","Date"),summarise,mean.rf.cm=mean(in.to.cm(Data.Value),na.rm=T))
rf.dat2$WY<-WY(rf.dat2$Date);#Water Year
rf.dat2$month<-as.numeric(format(rf.dat2$Date,"%m"))
rf.dat2$CY<-as.numeric(format(rf.dat2$Date,"%Y"))

# Monthly mean total rainfall for the period of record
rf.dat.mon<-ddply(rf.dat2,c("SITE","Region","CY","month","WY"),summarise,sum.cm=sum(mean.rf.cm,na.rm=T),N.val=N(mean.rf.cm))

#Ensure a full month is represented. 
rf.dat.mon.mean<-ddply(subset(rf.dat.mon,N.val>=28),c("Region","month"),summarise,mean.val=mean(sum.cm,na.rm=T),sd.val=sd(sum.cm,na.rm=T),N.val=N(sum.cm))
rf.dat.mon.mean<-merge(rf.dat.mon.mean,data.frame(month=c(5:12,1:4),month.order=1:12),"month")
rf.dat.mon.mean<-rf.dat.mon.mean[order(rf.dat.mon.mean$Region,rf.dat.mon.mean$month.order),]

# Annual total period of record
rf.dat.ann<-ddply(rf.dat2,c("SITE","Region","WY"),summarise,sum.cm=sum(mean.rf.cm,na.rm=T),N.val=N(mean.rf.cm))
rf.dat.ann<-subset(rf.dat.ann,N.val>364)
rf.dat.ann$sum.cm<-with(rf.dat.ann,ifelse(sum.cm==0,NA,sum.cm))

# Rainfall plots

cols<-wes_palette("Zissou1",6,"continuous")
regions.val=c("Orlando","Lake Okeechobee","West Palm","Southern WCA3A","ENP (SRS)", "ENP (Taylor Slough)")
# Monthly summary
xlim.val<-c(1,12);xmaj<-1:12;xmin<-1:12;xlab<-c(5:12,1:4)
ylim.val<-c(0,25);by.y<-5;ymaj<-seq(ylim.val[1],ylim.val[2],by.y);ymin<-seq(ylim.val[1],ylim.val[2],by.y/2)
plot(mean.val~month.order,rf.dat.mon.mean,ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA,type="n",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
polygon(c(1,6,6,1),c(-1,-1,30,30),border=NA,col=adjustcolor("skyblue",0.25))
text(3.5,24,"Wet Season",cex=0.75)
text(9.5,24,"Dry Season",cex=0.75)
for(i in 1:6){
  with(subset(rf.dat.mon.mean,Region==regions.val[i]),pt_line(month.order,mean.val,1,adjustcolor(cols[i],0.25),2,21,adjustcolor(cols[i],0.5),1.5,pt.lwd=0.1,pt.col=adjustcolor("black",0.5)))
}
axis_fun(1,line=-0.5,xmaj,xmin,xlab)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.25,"Monthly Mean Rainfall\n(cm month\u207B\u00B9)")
mtext(side=1,line=1.75,"Month")

# Annual summary
xlim.val<-c(1984,2019);by.x<-10;xmaj<-seq(xlim.val[1],xlim.val[2],by.x);xmin<-seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val<-c(0,300);by.y<-100;ymaj<-seq(ylim.val[1],ylim.val[2],by.y);ymin<-seq(ylim.val[1],ylim.val[2],by.y/2)
plot(sum.cm~WY,rf.dat.ann,ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:6){
  with(subset(rf.dat.ann,Region==regions.val[i]),pt_line(WY,sum.cm,2,adjustcolor(cols[i],0.25),1,21,adjustcolor(cols[i],0.5),1.25,pt.lwd=0.1,pt.col=adjustcolor("black",0.5)))
}
mean.rf=ddply(rf.dat.ann,"WY",summarise,mean.val=mean(sum.cm,na.rm=T),sd.val=sd(sum.cm,na.rm=T),N.val=N(sum.cm))
mean.rf$DOF=mean.rf$N.val-1
mean.rf$Tp=abs(qt(0.975,mean.rf$DOF))
mean.rf$U.CI=with(mean.rf,mean.val+(sd.val*Tp)/sqrt(N.val))
mean.rf$L.CI=with(mean.rf,mean.val-(sd.val*Tp)/sqrt(N.val))

with(mean.rf,shaded.range(WY,L.CI,U.CI,"grey",lty=1))
with(mean.rf,lines(WY,mean.val,col="black",lwd=1.5))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Annual Rainfall (cm yr\u207B\u00B9)")
mtext(side=1,line=1.75,"Water Year")


# End