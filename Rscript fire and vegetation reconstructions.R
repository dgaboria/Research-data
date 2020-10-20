# UTF-8
# Project : PhD in environmental sciences 2017-2020
# Author: Gaboriau Dorian - gaboriau.dorian@gmail.com - dorian.gaboriau@uqat.ca
# Thesis: Prediction of wildfire regimes up to 2100 in Northwest territories (NWT) and evaluation of the impacts for the Tlicho First Nation
# Chapter 2 : Reconstitution and caracterization of the past fire regime, vegetation and temperatures (Holocene period) in the central NWT, Canada 
# Last update 2020/07/20

rm(list = ls())  # Deleting variables from the environment R

##########################################################################################################
#Biomass Burning
##########################################################################################################
# Installation and importation of libraries
# Install paleofire package with github of Olivier Blarquez. #devtools::install_github("paleofire/paleofire")

library (devtools)
library(paleofire)
library(ggplot2)
library(GCD)

#Charge work environment
setwd("F:/2-En_cours/AXE_2_FEUX_PASSES/Charbons/5-REGFF-BB-FS/Code_FS/Code_FS_Copie_OK/")

# Data for EMILE LAKE - file corresponding to INPUTS of CharAnalysis

files=c("BBEmile.csv")
metadata=c("metadata1.csv")
mydata=pfAddData(files=files, metadata = metadata, type="CharAnalysis", sep = ";", dec = ".")

#Charcoal data transformation, background estimation and homogenization for unique to multiple series
TR1=pfTransform(add=mydata, method=c("MinMax","Box-Cox","Z-Score"))

#Produces a composite series from multiple charcoal records by using a robust locally weighted scatterplot smoother (LOWESS = locfit function from the locfit package and is applied repeatedly (nboot times) on bootstrapped charcoal sites samples. The records charcoal values are pre-binned prior to sites resampling (Daniau et al. (2012)).
COMP2=pfCompositeLF(TR1, hw=250, nboot=1000, tarAge=seq(-68,9530,1), conf = c(0.05, 0.95)) #hw = demi-largeur de fenetre

plot(COMP2,ylim=c(-2,2),main=c("Emile"))
EmileBB<-as.data.frame(COMP2$Result[,1:2])
colnames(EmileBB)=c("Age","LocFit")
#write.csv(EmileBB,file="C:/Users/client/Desktop/Emile.csv",row.names = FALSE)


# Data for IZAAC LAKE

files=c("BBIzaac.csv")
#colnames(files) = c("DepthTop", "DepthBottom", "AgeTop", "AgeBottom", "Volume", "Charcoal")
metadata=c("metadata2.csv")
mydata=pfAddData(files=files, metadata = metadata, type="CharAnalysis", sep = ";", dec = ".")
TR1=pfTransform(add=mydata, method=c("MinMax","Box-Cox","Z-Score"))
COMP3=pfCompositeLF(TR1, hw=250, nboot=1000, tarAge=seq(-68,9530,1), conf = c(0.05, 0.95))
plot(COMP3,ylim=c(-2,2),main=c("Izaac"))
IzaacBB<-as.data.frame(COMP3$Result[,1:2])
colnames(IzaacBB)=c("Age","LocFit")

#PARADIS LAKE

files=c("BBParadis.csv")
metadata=c("metadata3.csv")
mydata=pfAddData(files=files, metadata = metadata, type="CharAnalysis", sep = ";", dec = ".")
TR1=pfTransform(add=mydata, method=c("MinMax","Box-Cox","Z-Score"))
COMP5=pfCompositeLF(TR1, hw=250, nboot=1000, tarAge=seq(-68,9530,1), conf = c(0.05, 0.95))
plot(COMP5,ylim=c(-2,2),main=c("Paradis"))
ParadisBB<-as.data.frame(COMP5$Result[,1:2])
colnames(ParadisBB)=c("Age","LocFit")

#SAXON LAKE

files=c("BBSaxon.csv")
metadata=c("metadata4.csv")
mydata=pfAddData(files=files, metadata = metadata, type="CharAnalysis", sep = ";", dec = ".")
TR1=pfTransform(add=mydata, method=c("MinMax","Box-Cox","Z-Score"))
COMP4=pfCompositeLF(TR1, hw=250, nboot=1000, tarAge=seq(-68,9530,1), conf = c(0.05, 0.95))
plot(COMP4,ylim=c(-2,2),main=c("Saxon"))
SaxonBB<-as.data.frame(COMP4$Result[,1:2])
colnames(SaxonBB)=c("Age","LocFit")

#Merge BB of each lakes 
BB = merge(EmileBB, IzaacBB,  by='Age')
BB = merge(BB, ParadisBB, by = 'Age')
colnames(BB) = c("Age","Emile", "Izaac", "Paradis")
BB = merge(BB, SaxonBB, by = 'Age')
colnames(BB) = c("Age", "EmileBB", "IzaacBB", "ParadisBB", "SaxonBB")

#Measure the mean of biomass burning to have a regional value (BB_mean)
BB[,6]<-apply(BB[,2:ncol(BB)],1,mean,na.rm=T)
colnames(BB)=c("Age","EmileBB", "IzaacBB", "ParadisBB", "SaxonBB", "BBmean")
head(BB)

write.csv(BB,file="BB250.csv",row.names = FALSE)

#Plot BB of each lakes after having substracted the mean of reference period to each subsample from the file precedently created
dev.off()
#Use the file created
BB<-read.csv("F:/2-En_cours/AXE_2_FEUX_PASSES/Charbons/5-REGFF-BB-FS/Code_FS/Code_FS_Copie_OK/BB250ok.csv",h=T,sep=";")
colnames(BB) = c("Age", "EmileBB", "IzaacBB", "ParadisBB", "SaxonBB", "BBmean")
plot(BB$SaxonBB,xlim = rev(range(BB$Age)),type="l", ylim = c(-2,3), lwd = 2, ylab = "", col = "blue", cex.axis=1.3)
lines(BB$Age,BB$IzaacBB, col= "forestgreen", lwd = 2)
lines(BB$Age,BB$ParadisBB, col= "red", lwd = 2)
lines(BB$Age,BB$EmileBB, col= "black", lwd = 2)
abline(mean(BB$BBmean[69:569]),0) # = periode ca. 500 yr. BP to 0 (1950)

##########################################################################################################
#Fire frequency
##########################################################################################################
# Importation of libraries
library(devtools)
library(paleofire)

# Create a file for each lake, with the date of fire peaks (significant peaks)
Lakes_fires<-read.csv("F:/2-En_cours/AXE_2_FEUX_PASSES/Charbons/5-REGFF-BB-FS/Code_FS/Code_FS_Copie_OK/Lake_fire_frequency.csv",h=T,sep=";")
Emilefires<-Lakes_fires[,1]
Emilefires<-Emilefires[!is.na(Emilefires)]
Izaacfires<-Lakes_fires[,2]
Izaacfires<-Izaacfires[!is.na(Izaacfires)]
Paradisfires<-Lakes_fires[,3]
Paradisfires<-Paradisfires[!is.na(Paradisfires)]
Saxonfires<-Lakes_fires[,4]
Saxonfires<-Saxonfires[!is.na(Saxonfires)]
Saxonfires 

#Calcul du FF
fevent=c(Emilefires)

#Computes paleo-fire frequency for a set of fire events (or frequency from other events types, see examples) using a gaussian kernel density estimation procedure based on a defined bandwidth (see Mudelsee 2004 for details). Pseudo-replicated values are used to correct for edge bias, equivalent to "minimum slope" correction in Mann (2004).
ffEmile=kdffreq(fevent,up=-68,lo=9530, bandwidth = 500, nbboot=1000,alpha = 0.1, pseudo = FALSE)
#pseudo=FALSE is important because it disables the correction of Mann 2004 (pdf of paleofire package)

mat_ffEmile<-matrix(cbind(ffEmile$age,ffEmile$ff),ncol=2)
plot.kdffreq(ffEmile,ylim=c(0,0.015),xlim=c(10000,0),bty="n")
mat_ffEmile
colnames(mat_ffEmile)=c("Age","ff")

fevent=c(Izaacfires)
ffIzaac=kdffreq(fevent,up=-68,lo=9530, bandwidth = 500, nbboot=1000,alpha = 0.1, pseudo=FALSE)
mat_ffIzaac<-matrix(cbind(ffIzaac$age,ffIzaac$ff),ncol=2)
plot.kdffreq(ffIzaac,ylim=c(0,0.015),xlim=c(10000,0),bty="n")
mat_ffIzaac
colnames(mat_ffIzaac)=c("Age","ff")

fevent=c(Paradisfires)
ffParadis=kdffreq(fevent,up=-68,lo=9530,bandwidth = 500, nbboot=1000,alpha = 0.1,pseudo=FALSE)
mat_ffParadis<-matrix(cbind(ffParadis$age,ffParadis$ff),ncol=2)
plot.kdffreq(ffParadis,ylim=c(0,0.015),xlim=c(10000,0),bty="n")
mat_ffParadis
colnames(mat_ffParadis)=c("Age","ff")

fevent=c(Saxonfires)
ffSaxon=kdffreq(fevent,up=-68, lo=9530,bandwidth = 500, nbboot=1000,alpha = 0.1,pseudo=FALSE)
mat_ffSaxon<-matrix(cbind(ffSaxon$age,ffSaxon$ff),ncol=2)
plot.kdffreq(ffSaxon,ylim=c(0,0.015),xlim=c(10000,0),bty="n")
mat_ffSaxon
colnames(mat_ffSaxon)=c("Age","ff")

#Plot FF of each lakes 
plot(ffSaxon$ff, xlim = rev(range(ffSaxon$age)),type="l", ylab = "", xlab = "", ylim = c(0,0.010), lwd = 2, col = "blue", cex.axis=1.3)
lines(ffIzaac$age,ffIzaac$ff, col="forestgreen", lwd = 2)
lines(ffParadis$age,ffParadis$ff, col="red", lwd = 2)
lines(ffEmile$age,ffEmile$ff, col="black", lwd = 2)

#As dates are interpolated with the function kdffreq, they are not the same for each lake 
#Do a linear approximation all 10 years (or more or less)

mat_ffEmile<-approx(mat_ffEmile[,1],mat_ffEmile[,2],method="linear",xout=seq(-68,9530,by=1), rule = 2)
mat_ffIzaac<-approx(mat_ffIzaac[,1],mat_ffIzaac[,2],method="linear",xout=seq(-69,9530,by=1), rule = 2)
mat_ffParadis<-approx(mat_ffParadis[,1],mat_ffParadis[,2],method="linear",xout=seq(-68,9530,by=1), rule = 2)
mat_ffSaxon<-approx(mat_ffSaxon[,1],mat_ffSaxon[,2],method="linear",xout=seq(-68,9530,by=1), rule = 2)
FF<-as.data.frame(cbind(mat_ffEmile$x,mat_ffEmile$y, mat_ffIzaac$y, mat_ffParadis$y, mat_ffSaxon$y))
colnames(FF) = c("Age", "EmileFF", "IzaacFF", "ParadisFF", "SaxonFF")
dim(FF)
head(FF)
FF[,6]<-apply(FF[,2:5],1, mean)
head(FF)
colnames(FF)=c("Age","FFEmile","FFIzaac","FFParadis","FFSaxon", "FFmean")
abline(mean(FF$FFmean), 0)

write.csv(FF,file="FF500.csv",row.names = FALSE)

#Plot FF of each lakes after having substracted the mean of referecne period to each subsample from the file precedently created
dev.off()
#Use the file created
FF<-read.csv("F:/2-En_cours/AXE_2_FEUX_PASSES/Charbons/5-REGFF-BB-FS/Code_FS/Code_FS_Copie_OK/FF500ok.csv",h=T,sep=";")
colnames(FF) = c("Age", "EmileFF", "IzaacFF", "ParadisFF", "SaxonFF", "FFmean")
plot(FF$SaxonFF,xlim = rev(range(FF$Age)),type="l", ylim = c(-0.003,0.006), lwd = 2, ylab = "", col = "blue", cex.axis=1.3)
lines(FF$Age,FF$IzaacFF, col= "forestgreen", lwd = 2)
lines(FF$Age,FF$ParadisFF, col= "red", lwd = 2)
lines(FF$Age,FF$EmileFF, col= "black", lwd = 2)
abline(mean(FF$FFmean[69:569]),0)

#####################################################################################################################
#Measure RegFF
#####################################################################################################################
dev.off()
par(mfrow = c(4,1), mar = c(0.1,7,4,0.1))

#Import Libraries
library(boot)

#File with fire frequency for each lake and mean for each year (corrected with the difference with the mean of ca. 3000-0 yr. BP)
FFSap<-read.csv("FF500ok.csv",h=T,sep=";")
FFSap<-na.omit(FFSap)
FFSap_plot<-FFSap[,c("FFEmile","FFIzaac","FFParadis","FFSaxon")]
my.data<-as.data.frame(t(FFSap_plot))
R<-999
l=length(my.data[1,])
bootCI<-as.data.frame(matrix(ncol=3,nrow=l))
boot.mean<-function(y,m)
{z<-mean(y[m])
z}
for(i in 1:l) {
  y<-my.data[,i]
  y<-na.omit(y)
  boot.obj<-boot(y,boot.mean,R)
  IC.norm<-boot.ci(boot.obj,conf=0.90,type="norm")
  bootCI[i,]<-IC.norm$normal
}
FFBootCISap<-cbind(FFSap$Age,FFSap$FFmean,bootCI)
colnames(FFBootCISap)=c("Age","Mean","CI","Infe","Supe")
FFmean_smoothSap = smooth.spline(FFBootCISap$Age,FFBootCISap$Mean, spar=0.15)
FFmean_infSap = smooth.spline(FFBootCISap$Age, FFBootCISap$Infe, spar=0.15)
FFmean_supSap = smooth.spline(FFBootCISap$Age, FFBootCISap$Supe, spar=0.15)
plot(FFmean_infSap,type="l",xlim=c(9530,-68),ylim=c(-0.004,0.006), axes=F, xlab = " ", xaxt="n", lty = 1, lwd = 2, font.axis=3, cex.axis = 1, ylab = " ", col=c("yellow"), bty="n")
lines(FFmean_supSap,col=c("yellow"), lwd = 2)
polygon(c(FFmean_smoothSap$x,rev(FFmean_smoothSap$x)),c(FFmean_supSap$y,rev(FFmean_infSap$y)),col="yellow", border = NA)
linemeanFF<-mean(FFmean_smoothSap$y[69:569])
lines(FF$Age,FF$IzaacFF, col="orange", lwd = 2.5, lty = 2)
lines(FF$Age,FF$ParadisFF, col="brown", lwd = 2.5, lty = 2)
lines(FF$Age,FF$EmileFF, col="forestgreen", lwd = 2.5, lty = 2)
lines(FF$Age,FF$SaxonFF, col="blue", lwd = 2.5, lty = 2)
lines(FF$Age,FF$FFmean, col= "red", lwd = 3)
abline(a=NULL,b=NULL,h=linemeanFF,v=NULL, col = "black", lwd = 1, lty = 1)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)
axis(3,at=c(-68,0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5)

#The next two files contain the average of RegBB and RegFF and the INF and SUP bounds of the Confidence Intervals.
write.csv(BBBootCISap,"BBBootCISap.csv",row.names=F)
write.csv(FFBootCISap,"FFBootCISap.csv",row.names=F)

#####################################################################################################################
#Measure RegBB
#####################################################################################################################

#Import libraries
library(boot)

#Set environment
setwd("F:/2-En_cours/AXE_2_FEUX_PASSES/Charbons/5-REGFF-BB-FS/Code_FS/Code_FS_Copie_OK/") 

#File with Burned Biomass for each lake and mean for each year (corrected with the difference with the mean of ca. 3000-0 yr. BP)
BBSap<-read.csv("BB250ok.csv",h=T,sep=";")

#Keep just lakes without mean for bootstrap
BBSap_plot<-BBSap[,c("EmileBB","IzaacBB","ParadisBB","SaxonBB")]

#Confidence interval around RegBB with Boostrap
my.data<-as.data.frame(t(BBSap_plot)) # you return your table (format required to use the boot function)
R<-999  # 1000 is the number of resamples you want in your bootstrap
l=length(my.data[1,])
bootCI<-as.data.frame(matrix(ncol=3,nrow=l))  #You create your empty matrix to accommodate the output data
boot.mean<-function(y,m)       #You create your function that will calculate your bootstrap around your average
{z<-mean(y[m])
z}

for(i in 1:l) {
  y<-my.data[,i]
  y<-na.omit(y)
  boot.obj<-boot(y,boot.mean,R)
  IC.norm<-boot.ci(boot.obj,conf=0.90,type="norm")     #you can change confidence interval (here conf=0.90)
  bootCI[i,]<-IC.norm$normal
}

#To add the IC to the dataframe

colnames(BBSap) = c("Age", "EmileBB", "IzaacBB", "ParadisBB", "SaxonBB", "BBmean")
BBBootCISap<-cbind(BBSap$Age,BBSap$BBmean,bootCI)
colnames(BBBootCISap)=c("Age","Mean","CI","Infe","Supe")
BBmean_smoothSap = smooth.spline(BBBootCISap$Age,BBBootCISap$Mean, spar=0.15)   #smoothing of the data

#Plot RegBB and IC

BBmean_infSap = smooth.spline(BBBootCISap$Age, BBBootCISap$Infe, spar=0.15)
BBmean_supSap = smooth.spline(BBBootCISap$Age, BBBootCISap$Supe, spar=0.15)
plot(BBmean_infSap,type="l",xlim=c(9530,-68),ylim=c(-1.45,3),col=c("yellow"), axes=F, font.axis=2, xlab = " ", ylab = " ", cex.axis=0.8, lty=1, bty="n", lwd = 2, xaxt="n")
lines(BBmean_supSap,col=c("yellow"), lwd = 2, lty = 1)
polygon(c(BBmean_smoothSap$x,rev(BBmean_smoothSap$x)),c(BBmean_supSap$y,rev(BBmean_infSap$y)),col="yellow", border = NA)
lines(BB$Age,BB$IzaacBB, col= "orange", lwd = 2.5, lty =2)
lines(BB$Age,BB$ParadisBB, col= "brown", lwd = 2.5, lty = 2)
lines(BB$Age,BB$EmileBB, col= "forestgreen", lwd = 2.5, lty = 2)
lines(BB$Age,BB$SaxonBB, col= "blue", lwd = 2.5, lty = 2)
linemeanBB<-mean(BBmean_smoothSap$y[69:569])
lines(BB$Age,BB$BBmean, col= "red", lwd = 3)
abline(a=NULL,b=NULL,h=linemeanBB,v=NULL, col = "black", lwd = 1, lty = 1)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)

##################################################################################################################
#### FS index
##################################################################################################################

BB<-read.csv("BBBootCISap.csv",h=T,sep=",")
FF<-read.csv("FFBootCISap.csv",h=T,sep=",")
Tot<-merge(BB,FF,by="Age")
BB<-cbind(Tot[,1],Tot[,2],Tot[,3],Tot[,4],Tot[,5])
colnames(BB)<-c("Age","Mean","CI","Infe","Supe")
FF<-cbind(Tot[,1],Tot[,6],Tot[,7],Tot[,8],Tot[,9])
colnames(FF)<-c("Age","Mean","CI","Infe","Supe")

RegBB<-data.frame(BB[,1],BB[,2])
colnames(RegBB)<-c("Age","BBmean")
RegBBmax<-max(RegBB$BBmean)
RegBBmin<-min(RegBB$BBmean)

#Here is a calcul to center your RegBB around the average.# It's not important if you calculate a single RegBB (it won't change your trend but only your RegBB value gradient) but it's very important if you make a comparison of several RegBB/RegFF and/or FS. If this is the case, you have to center all your RegBBs on the average of all your RegBBs because otherwise you will compare results that are not really comparable. 
RegBB[,3]<-c((RegBB[,2]-RegBBmin)/(RegBBmax-RegBBmin))
colnames(RegBB)<-c("Age","BBmean","RegBBrescale")

#We add a +1 to have only data higher than 1 (useful for the continuation of the calculation and mention in Olivier's method in the article of the PNAS)
RegBB[,4]<-RegBB$RegBBrescale+1
colnames(RegBB)<-c("Age","BBmean","RegBBrescale","RegBBb")

#Same for RegFF
RegFF<-data.frame(FF[,1],FF[,2])
colnames(RegFF)<-c("Age","FFmean")                
RegFFmax<-max(na.omit(RegFF$FFmean))
RegFFmin<-min(na.omit(RegFF$FFmean))
RegFF[,3]<-c((RegFF[,2]-RegFFmin)/(RegFFmax-RegFFmin))
colnames(RegFF)<-c("Age","FFmean","RegFFrescale")
RegFF[,4]<-RegFF$RegFFrescale+1
colnames(RegFF)<-c("Age","FFmean","RegFFrescale","RegFFb")

#Calcul of FS index
FS<-merge(RegBB,RegFF,by="Age",all=T)
FS[,8]<-(FS$RegBBb/FS$RegFFb)  
FSspline<-FS[1:nrow(FS),c(1,8)]
colnames(FSspline)<-c("Age","FSindex")
write.csv(FSspline,file="FS500.csv",row.names = FALSE)

#Import file corrected with difference with the reference period ca. 3000-0 yr. BP
FSspline = read.csv("FS500ok.csv",h=T,sep=";")

#Then you can make a smooth to make a graph. You can use other methods than smooth.spline. 
FSspline<-na.omit(FSspline)
smoothingFSSpline = smooth.spline(FSspline$Age, FSspline$FSindex, spar=0.15)

MAXmean<-RegBBmax
MINmean<-RegBBmin

#Calcul IC around FS index, this time with the confidence intervals around your RegBB and RegFF
BBBootCISap<-read.csv("BBBootCISap.csv",h=T,sep=",")
FFBootCISap<-read.csv("FFBootCISap.csv",h=T,sep=",")

FS_CI<-merge(BBBootCISap,FFBootCISap,by="Age")
FS_CI<-cbind(FS_CI[,1],FS_CI[,2],FS_CI[,4],FS_CI[,5],FS_CI[,6],FS_CI[,8],FS_CI[,9])
colnames(FS_CI)<-c("Age","BBmean","BBinf","BBsup","FFmean","FFinf","FFsup")
write.csv(FS_CI,"FS_CI.csv",row.names=F)

RegBB<-c()
RegFF<-c()
RegBB<-cbind(FS_CI[,1],FS_CI[,3],FS_CI[,2])
RegBB<-data.frame(RegBB)
colnames(RegBB)<-c("Age","BBinf","BBmean")
RegBBmax<-max(RegBB$BBmean)
RegBBmin<-min(RegBB$BBmean)
RegBB[,4]<-c((RegBB[,2]-RegBBmin)/(RegBBmax-RegBBmin))
colnames(RegBB)<-c("Age","BBinf","BBmean","RegBBrescale")
RegBB[,5]<-RegBB$RegBBrescale+1
colnames(RegBB)<-c("Age","BBinf","BBmean","RegBBrescale","RegBBb")

RegFF<-cbind(FS_CI[,1],FS_CI[,7],FS_CI[,5])
RegFF<-data.frame(RegFF)
colnames(RegFF)<-c("Age","FFsup","FFmean")                
RegFFmax<-max(RegFF$FFmean)
RegFFmin<-min(RegFF$FFmean)
RegFF[,4]<-c((RegFF[,2]-RegFFmin)/(RegFFmax-RegFFmin))
colnames(RegFF)<-c("Age","FFsup","FFmean","RegFFrescale")
RegFF[,5]<-RegFF$RegFFrescale+1
colnames(RegFF)<-c("Age","FFsup","FFmean","RegFFrescale","RegFFb")

FS_inf<-merge(RegBB,RegFF,by="Age",all=T)
FS_inf<-(FS_inf$RegBBb/FS_inf$RegFFb)

RegBB<-c()
RegFF<-c()
RegBB<-cbind(FS_CI[,1],FS_CI[,4],FS_CI[,2])
RegBB<-data.frame(RegBB)
colnames(RegBB)<-c("Age","BBsup","BBmean")
RegBBmax<-max(RegBB$BBmean)
RegBBmin<-min(RegBB$BBmean)
RegBB[,4]<-c((RegBB[,2]-RegBBmin)/(RegBBmax-RegBBmin))
colnames(RegBB)<-c("Age","BBsup","BBmean","RegBBrescale")
RegBB[,5]<-RegBB$RegBBrescale+1
colnames(RegBB)<-c("Age","BBsup","BBmean","RegBBrescale","RegBBb")

RegFF<-cbind(FS_CI[,1],FS_CI[,6],FS_CI[,5])
RegFF<-data.frame(RegFF)
colnames(RegFF)<-c("Age","FFinf","FFmean")                
RegFFmax<-max(RegFF$FFmean)
RegFFmin<-min(RegFF$FFmean)
RegFF[,4]<-c((RegFF[,2]-RegFFmin)/(RegFFmax-RegFFmin))
colnames(RegFF)<-c("Age","FFinf","FFmean","RegFFrescale")
RegFF[,5]<-RegFF$RegFFrescale+1
colnames(RegFF)<-c("Age","FFinf","FFmean","RegFFrescale","RegFFb")

FS_sup<-merge(RegBB,RegFF,by="Age",all=T)
FS_sup<-(FS_sup$RegBBb/FS_sup$RegFFb)

FS_CI<-cbind(FS_CI[,1],FS_inf,FS_sup)
colnames(FS_CI)<-c("Age","Infe","Supe")

write.csv(FS_CI,"FS_CI_2.csv",row.names=F)

#Import file corrected with differences between each subsamples and the mean of the reference period ca. 3000-0 yr. BP
FS_CI<-read.csv("FS_CI_2ok.csv",h=T,sep=";")
FS<-read.csv("FS500ok.csv",h=T,sep=";")
FS<-merge(FS,FS_CI,by="Age") 
FS<-as.data.frame(FS)

smoothingFSSap = smooth.spline(FSspline[,1], FSspline[,2] , spar=0.15) 
plot(smoothingFSSap,type="l",xlim=c(9530,-68),ylim=c(-0.5,1.5),col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=2,cex.axis=1.3,bty="n", xaxt="n")
FS_infSap = smooth.spline(FS[,1], FS[,3], spar=0.15)
FS_supSap = smooth.spline(FS[,1], FS[,4], spar=0.15)
lines(FS_infSap,col=c("yellow"), lwd = 2)
lines(FS_supSap,col=c("yellow"), lwd = 2)
polygon(c(smoothingFSSap$x, rev(smoothingFSSap$x)), c(FS_supSap$y, rev(FS_infSap$y)),col = "yellow", border = NA)
lines(smoothingFSSap, col=c("red"), lwd = 3)
linemeanFS<-mean(smoothingFSSap$y[69:569])
abline(a=NULL,b=NULL,h=linemeanFS,v=NULL, col = "black", lwd = 1, lty = 1)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)

# Data Reg BB, RegFF, FSindex. Merge in a dataframe

firemetrics = as.data.frame(cbind(BBBootCISap$Age, BBBootCISap$Mean, FFBootCISap$Mean, FSspline$FSindex))
colnames(firemetrics) = c("Age","RegBB", "RegFF", "FSindex")
write.csv(firemetrics,file="F:/2-En_cours/AXE_2_FEUX_PASSES/firemetrics.csv",row.names = FALSE)

################################################################
################################################################
#VEGETATION reconstructions
################################################################
################################################################

# Graphic for all taxa POLLEN INFLUX DIAGRAM

setwd("F://2-En_cours/AXE_2_FEUX_PASSES/Charbons/3-POLLENS/data_diagram/")

# Import data
ma.pollen <- read.csv("Pollen_influx.csv", header=TRUE, sep=";", check.names = FALSE)
dev.off()

#PLOT TREES & SHRUBS
ma.pollen_1 = as.data.frame(cbind(as.numeric(ma.pollen$`ï»¿age`), ma.pollen$PAR, ma.pollen$Picea, ma.pollen$Betula, ma.pollen$Pinus, ma.pollen$Populus, ma.pollen$Alnus.crispa, ma.pollen$Alnus.rugosa, ma.pollen$Juniperus, ma.pollen$Salix, ma.pollen$Larix))
colnames(ma.pollen_1) = c("Age", "PAR", "Picea", "Betula", "Pinus", "Populus", "Alnus crispa", "Alnus rugosa", "Juniperus","Salix", "Larix")
ma.pollen_1 = ma.pollen_1[1:148,]

# Define colour scheme
p.col <- c(rep("red", times = 1), rep("darkred", times = 4), rep("forestgreen", times = 5))
pol.plot1 = strat.plot(ma.pollen_1[,2:11],exag = T, exag.mult = 5, exag.alpha=0.2, yvar = ma.pollen_1$Age, col.poly.line="black", cex.xlabel = 1.2, scale.minmax = T, xSpace=0.015, srt.xlabel = 70, scale.percent=TRUE, y.rev = T, plot.line=T, col.line = p.col, plot.poly=TRUE, plot.bar=F, col.poly=p.col)

#PLOT PLANTS & AQUATICS
ma.pollen_2 = as.data.frame(cbind(as.numeric(ma.pollen$`ï»¿age`), ma.pollen$Cyperaceae, ma.pollen$Ericaceae, ma.pollen$Myrica, ma.pollen$Lycopodium, ma.pollen$Artemisia, ma.pollen$Poaceae, ma.pollen$Pediastrum, ma.pollen$Nuphar, ma.pollen$Potamogeton))
colnames(ma.pollen_2) = c("Age", "Cyperaceae", "Ericaceae","Myrica", "Lycopodium","Artemisia", "Poaceae", "Pediastrum", "Nuphar", "Potamogeton")
ma.pollen_2 = ma.pollen_2[1:148,]

# Define colour scheme
p.col <- c(rep("gold2", times = 6), rep("blue", times = 3))
pol.plot2 = strat.plot(ma.pollen_2[,2:10], exag = T, exag.mult = 5, exag.alpha=0.2, yvar = ma.pollen_1$Age, col.poly.line="black", cex.xlabel = 1.2, scale.minmax = T, xSpace=0.015, srt.xlabel = 70, scale.percent=TRUE, y.rev = T, plot.line=T, col.line = p.col, plot.poly=TRUE, plot.bar=F, col.poly=p.col)

# POLLEN DIAGRAM (%)
library(rioja)

setwd("F:/2-En_cours/AXE_2_FEUX_PASSES/Charbons/3-POLLENS/data_diagram/")

# Import data
ma.pollen <- read.csv("pollen_percentages.csv", header=TRUE, sep=";", check.names = FALSE)
txsedim <- read.csv("Sedimentation_rate.csv", header=TRUE, sep=";", check.names = FALSE)

dev.off()

x  = pol.plot <- strat.plot(ma.pollen[,4:23], scale.minmax = T, x.pc.inc = 10, yvar=ma.pollen$`ï»¿Age`, xLeft=0.15, exag.mult=5, y.rev=TRUE, exag = T, col.exag = "auto", exag.alpha=0.2, plot.line=TRUE, plot.poly=T, plot.bar=F, col.poly = "black", lwd.bar=10, sep.bar=TRUE, scale.percent=TRUE, xSpace=0.01, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2)
x  = pol.plot <- strat.plot(ma.pollen[,4:23], scale.minmax = T, x.pc.inc = 10, yvar=ma.pollen$`Depth (cm)`,y.tks = , xLeft=0.15, exag.mult=5, y.rev=TRUE, exag = T, col.exag = "auto", exag.alpha=0.2, plot.line=TRUE, plot.poly=T, plot.bar=F, col.poly = "black", lwd.bar=10, sep.bar=TRUE, scale.percent=TRUE, xSpace=0.01, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2)

plot(txsedim$`Depth (cm)`, txsedim$nbr_yr, type = "l")

library(ggplot2)
ggplot(txsedim, aes(seq(length=nrow(txsedim)), txsedim$nbr_yr)) + geom_path() + #Ploting
  scale_y_continuous(name= "Number of failures") +
  scale_x_continuous(name= "Operations performed")

##########################################################################################################################
#Temperatures
##########################################################################################################################

par(mfrow = c(5,1), mar = c(0.3,4,4,1))

##### UPITER ET AL., 2014 (chironomid inferred July temperature) - ncdc.noaa.gov - 2 transfert functions (Barley et al., 2006 & Porinchu et al., 2009)

FStemp<-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/Data/Data_Jesse_Vermaire_Carleton_NWT/Data_transfer_Porinchu_2009.csv",h=T,sep=";")
head(FStemp)
colnames(FStemp) = c("Age", "Temp", "sup", "inf")

smoothingFStemp = smooth.spline(FStemp$Age, FStemp$Temp, spar=0.15) 
plot(smoothingFStemp,type="l",xlim=c(10000,-68),ylim=c(7.5,14.5),col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=2,cex.axis=1.3,bty="n", xaxt="n")
FStemp_infSap = smooth.spline(FStemp$Age, FStemp$inf, spar=0.15)
FStemp_supSap = smooth.spline(FStemp$Age, FStemp$sup, spar=0.15)
polygon(c(smoothingFStemp$x, rev(smoothingFStemp$x)), c(FStemp$sup, rev(FStemp$inf)),
        col = "aquamarine4", border = NA)
lines(smoothingFStemp, col=c("black"), lwd = 2)
linemeanFStemp<-mean(FStemp$Temp)
abline(a=NULL,b=NULL,h=linemeanFStemp,v=NULL, col = "red", lwd = 2, lty = 2)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)
axis(3,at=seq(0,10000,1000), labels=c("0","1","2","3","4","5","6", "7", "8", "9", "10") , lty = 1, lwd = 2, font.axis=2, cex.axis = 2)


##### TEMPERATURE FROM PORTER ET AL., 2019

FStemp2<-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Data_climat/Data/Data_Porter_Yukon/Temp_Porter.csv",h=T,sep=";")
head(FStemp2)
colnames(FStemp2) = c("Age", "Temp", "sup", "inf")
FStemp2 = subset(FStemp2, FStemp2$Age < 10000)

smoothingFStemp2 = smooth.spline(FStemp2$Age, FStemp2$Temp, spar=0.15) 
plot(smoothingFStemp2,type="l",xlim=c(10000,-68),ylim=c(-2,3),col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=2,cex.axis=1.3,bty="n", xaxt="n")
FStemp_infSap2 = smooth.spline(FStemp2$Age, FStemp2$inf, spar=0.15)
FStemp_supSap2 = smooth.spline(FStemp2$Age, FStemp2$sup, spar=0.15)
polygon(c(smoothingFStemp2$x, rev(smoothingFStemp2$x)), c(FStemp2$sup, rev(FStemp2$inf)),
        col = "aquamarine4", border = NA)
lines(smoothingFStemp2, col=c("black"), lwd = 2)
abline(a=NULL,b=NULL,h=0,v=NULL, col = "red", lwd = 2, lty = 2)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)

# #Approximate temperature 2
# approxtemp2 <-approx(FStemp2$Age,FStemp2$Temp ,method="linear",xout=seq(-68,9760,by=1), rule = 2)

##### LECAVALIER ET AL., 2017 - Arctic air temperature reconstruction 

FStemp3<-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Data_climat/Data/Data_lecavalier/Temp3_lecavalier.csv",h=T,sep=";")
head(FStemp3)
colnames(FStemp3) = c("Age", "Temp", "sup", "inf")
head(FStemp3)
FStemp3 = subset(FStemp3, FStemp3$Age < 10000)

smoothingFStemp3 = smooth.spline(FStemp3$Age, FStemp3$Temp, spar=0.15) 
plot(smoothingFStemp3,type="l",xlim=c(10000,-68),ylim=c(-2,8),col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=2,cex.axis=1.3,bty="n", xaxt="n")
FStemp_infSap3 = smooth.spline(FStemp3$Age, FStemp3$inf, spar=0.15)
FStemp_supSap3 = smooth.spline(FStemp3$Age, FStemp3$sup, spar=0.15)
polygon(c(smoothingFStemp3$x, rev(smoothingFStemp3$x)), c(FStemp3$sup, rev(FStemp3$inf)),
        col = "aquamarine4", border = NA)
lines(smoothingFStemp3, col=c("black"), lwd = 2)
linemeanFStemp3<-mean(FStemp3$Temp)
abline(a=NULL,b=NULL,h=linemeanFStemp3,v=NULL, col = "red", lwd = 2, lty = 2)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)

##### KOBASHI 2017

FStemp4<-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Data_climat/Data/Data_Kobashi/Temp4_kobashi2017.csv",h=T,sep=";")
head(FStemp4)
colnames(FStemp4) = c("Age", "Temp", "sup", "inf")
FStemp4 = subset(FStemp4, FStemp4$Age < 10000)
head(FStemp4)

smoothingFStemp4 = smooth.spline(FStemp4$Age, FStemp4$Temp, spar=0.15) 
plot(smoothingFStemp4,type="l",xlim=c(10000,-68),ylim=c(-36,-25),col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=2,cex.axis=1.3,bty="n", xaxt="n")
FStemp_infSap4 = smooth.spline(FStemp4$Age, FStemp4$inf, spar=0.15)
FStemp_supSap4 = smooth.spline(FStemp4$Age, FStemp4$sup, spar=0.15)
polygon(c(smoothingFStemp4$x, rev(smoothingFStemp4$x)), c(FStemp4$sup, rev(FStemp4$inf)),
        col = "aquamarine4", border = NA)
lines(smoothingFStemp4, col=c("black"), lwd = 2)
linemeanFStemp4<-mean(FStemp4$Temp)
abline(a=NULL,b=NULL,h=linemeanFStemp4,v=NULL, col = "red", lwd = 2, lty = 2)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)


##### Ajout temperatures V ##### VIAU 2006
FStemp5<-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Data_climat/Data/Data_Viau/Viau.csv",h=T,sep=";")
head(FStemp5)

colnames(FStemp5) = c("Age", "Temp")
FStemp5 = subset(FStemp5, FStemp5$Age < 10000)
head(FStemp5)

smoothingFStemp5 = smooth.spline(FStemp5$Age, FStemp5$Temp, spar=0.15) 
plot(smoothingFStemp5,type="l",xlim=c(10000,-68),ylim=c(12,16),col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=2,cex.axis=1.5,bty="n", xaxt="n")
linemeanFStemp5<-mean(FStemp5$Temp)
abline(a=NULL,b=NULL,h=linemeanFStemp5,v=NULL, col = "red", lwd = 2, lty = 2)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)

##SCALING
par(mfrow = c(5,1), mar = c(0.2,2,2,0.7), oma = c(1, 0.1, 0, 0.5))

FStemp1 = FStemp[,1:2]
approxtemp1 <-approx(FStemp1$Age,FStemp1$Temp ,method="linear",xout=seq(-68,6000,by=1), rule = 2)
approxtemp1$y = scale(approxtemp1$y, scale = T, center = T)
plot(approxtemp1$y, type = "l")
mean(approxtemp1$y)
sd(approxtemp1$y)

FStemp2 = FStemp2[,1:2]
approxtemp2 <-approx(FStemp2$Age,FStemp2$Temp ,method="linear",xout=seq(-68,10000,by=1), rule = 2)
approxtemp2$y = scale(approxtemp2$y, scale = T, center = T)
plot(approxtemp2$y, type = "l")
mean(approxtemp2$y)
sd(approxtemp2$y)

FStemp3 = FStemp3[,1:2]
approxtemp3 <-approx(FStemp3$Age,FStemp3$Temp ,method="linear",xout=seq(-68,10000,by=1), rule = 2)
approxtemp3$y = scale(approxtemp3$y, scale = T, center = T)
plot(approxtemp3$y, type = "l")
mean(approxtemp3$y)
sd(approxtemp3$y)

FStemp4 = FStemp4[,1:2]
approxtemp4 <-approx(FStemp4$Age,FStemp4$Temp ,method="linear",xout=seq(-68,10000,by=1), rule = 2)
approxtemp4$y = scale(approxtemp4$y, scale = T, center = T)
plot(approxtemp4$y, type = "l")
mean(approxtemp4$y)
sd(approxtemp4$y)

FStemp5 = FStemp5[,1:2]
approxtemp5 <-approx(FStemp5$Age,FStemp5$Temp ,method="linear",xout=seq(-68,10000,by=1), rule = 2)
approxtemp5$y = scale(approxtemp5$y, scale = T, center = T)
plot(approxtemp5$y, type = "l")
mean(approxtemp5$y)
sd(approxtemp5$y)

#write.csv(approxtemp1,file="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/TEMP1.csv",row.names = FALSE)
#write.csv(approxtemp2,file="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/TEMP2.csv",row.names = FALSE)
#write.csv(approxtemp3,file="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/TEMP3.csv",row.names = FALSE)
#write.csv(approxtemp4,file="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/TEMP4.csv",row.names = FALSE)
#write.csv(approxtemp5,file="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/TEMP5.csv",row.names = FALSE)

#Import pooled data (scaled)

dev.off()
#File with substraction between each subsample and mean of the reference period

TEMP <- read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/Mean_temp_scaled_0-500ka.csv", row.names = 1, header=TRUE, sep=";", check.names=F)
#write.csv(TEMP,file="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/Mean_temp_scaled2.csv",row.names = FALSE)

plot(rev(TEMP$MEAN),type="l",col=c("black"),ylab = " ", axes=F, font.axis=2, lwd=1,cex.axis=1.3,bty="n", xaxt="n")
linemean<-mean(TEMP$MEAN[69:569])
abline(a=NULL,b=NULL,h=linemean,v=NULL, col = "red", lwd = 2, lty = 2)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.5, las = 1)
axis(3,at=seq(0,10000,1000), labels=c("10","9","8","7","6","5","4", "3","2", "1", "0") , lty = 1, lwd = 2, font.axis=2, cex.axis = 2)

#Bootstrap around scaled temperature
library(boot)
dev.off()
head(TEMP)

TempSap<-TEMP[,1:5]
tempSap_plot<-TempSap[,c("TEMP1","TEMP2","TEMP3","TEMP4", "TEMP5")]  
my.data<-as.data.frame(t(tempSap_plot))
R<-999
l=length(my.data[1,])
bootCI<-as.data.frame(matrix(ncol=3,nrow=l))
boot.mean<-function(y,m)  
{z<-mean(y[m])
z}

for(i in 1:l) {
  y<-my.data[,i]
  y<-na.omit(y)
  boot.obj<-boot(y,boot.mean,R)
  IC.norm<-boot.ci(boot.obj,conf=0.90,type="norm")
  bootCI[i,]<-IC.norm$normal
}

#Add Interval Confidence in first dataframe
tempBootCISap<-cbind(as.numeric(rownames(TempSap)), TEMP$MEAN, bootCI)
colnames(tempBootCISap)=c("Age","Mean","CI","Infe","Supe")

moyenne = mean((tempBootCISap$Mean[69:569]))
summary(tempBootCISap)

Tempmean_smoothSap = tempBootCISap[,1:2]
Tempmean_smoothSap = smooth.spline(tempBootCISap$Mean, spar=0.15)   #Tu smooth tes données

#Plot TEMP and Interval of Confidence
par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Tempmean_smoothSap$x, Tempmean_smoothSap$y,type="l",xlim=c(10000,-68),ylim=c(-2,4),col=c("black"), axes=F, font.axis=2, xlab = " ", ylab = " ", cex.axis=0.8, lty=1, bty="n", lwd = 2, xaxt="n")
Tempmean_infSap = smooth.spline(tempBootCISap$Age, tempBootCISap$Infe, spar=0.15)
Tempmean_supSap = smooth.spline(tempBootCISap$Age, tempBootCISap$Supe, spar=0.15)
polygon(c(Tempmean_smoothSap$x,rev(Tempmean_smoothSap$x)),c(Tempmean_supSap$y,rev(Tempmean_infSap$y)),col="yellow2", border = NA)
lines(Tempmean_smoothSap, col=c("red"), lwd = 3)
linemeanTemp<-moyenne
abline(a=NULL,b=NULL,h=linemeanTemp,v=NULL, col = "black", lwd = 1, lty = 1)
axis(2, at=NULL, labels=TRUE, lty = 1, lwd = 2, font.axis=2, cex.axis = 1.3, las = 1, col = "black")
axis(3,at=c(-68,0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), lty = 1, lwd = 2, font.axis=2, cex.axis = 1.3)

#PROCESS CORRELATIONS BINCOR

rm(list = ls())  # Deleting variables from the environment R
dev.off() 

#Data_feu
firemetrics<-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/firemetrics.csv",h=T,sep=",")

#Data climat
climat <-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/Mean_temp_scaled_0-500ka.csv",h=T,sep=";")
#pp = read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/DATA_climat/pp_viau_2009.csv",h=T,sep=";")
#Data pollen (PAR)

pollen <-read.csv("E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Charbons/3-POLLENS/data_diagram/Pollen_influx.csv",h=T,sep=";")

# Correlations
#install.packages("BINCOR")
library(BINCOR)
#install.packages("pracma")
library(pracma)
###TEMP et fire metrics

RegFF = cbind(firemetrics$Age[1:9599], firemetrics$RegFF[1:9599])
RegBB = cbind(firemetrics$Age[1:9599], firemetrics$RegBB[1:9599])
FS = cbind(firemetrics$Age[1:9599], firemetrics$FSindex[1:9599])
TEMP = cbind(climat$ï..AGE, climat$MEAN)
TEMP8000 = cbind(climat$ï..AGE[1:8069], climat$MEAN[1:8069])
PAR = cbind(pollen$ï..age, pollen$PAR)
Picea = cbind(pollen$ï..age, pollen$Picea)
Betula = cbind(pollen$ï..age, pollen$Betula)
Alnus.c = cbind(pollen$ï..age, pollen$Alnus.crispa)
Alnus.r = cbind(pollen$ï..age, pollen$Alnus.rugosa)
Pinus = cbind(pollen$ï..age, pollen$Pinus)
Juniperus = cbind(pollen$ï..age, pollen$Juniperus)
Populus = cbind(pollen$ï..age, pollen$Populus)
Myrica = cbind(pollen$ï..age, pollen$Myrica)
Poaceae = cbind(pollen$ï..age, pollen$Poaceae)
Salix = cbind(pollen$ï..age, pollen$Salix)
Larix = cbind(pollen$ï..age, pollen$Larix)
Lycopodium = cbind(pollen$ï..age, pollen$Lycopodium)
Artemisia = cbind(pollen$ï..age, pollen$Artemisia)
Cyperaceae = cbind(pollen$ï..age, pollen$Cyperaceae)
Ericaceae = cbind(pollen$ï..age, pollen$Ericaceae)
Nuphar = cbind(pollen$ï..age, pollen$Nuphar)
Pediastrum = cbind(pollen$ï..age, pollen$Pediastrum)
Potamogeton = cbind(pollen$ï..age, pollen$Potamogeton)

# REGFF
test1 = bin_cor(RegFF, TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF")
binnedts<- test1$Binned_time_series
bin_ts1 <- na.omit(test1$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test1$Binned_time_series[,c(1,3)])
plot_ts(RegFF, TEMP, bin_ts1, bin_ts2, "RegFF", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test1 = bin_cor(RegFF, TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF")
binnedts<- test1$Binned_time_series
bin_ts1 <- na.omit(test1$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test1$Binned_time_series[,c(1,3)])
plot_ts(RegFF, TEMP, bin_ts1, bin_ts2, "RegFF", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test1 = bin_cor(RegFF, TEMP8000, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF")
binnedts<- test1$Binned_time_series
bin_ts1 <- na.omit(test1$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test1$Binned_time_series[,c(1,3)])
plot_ts(RegFF, TEMP8000, bin_ts1, bin_ts2, "RegFF", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test1 = bin_cor(RegFF, PP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_pp_RegFF")
# binnedts<- test1$Binned_time_series
# bin_ts1 <- na.omit(test1$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test1$Binned_time_series[,c(1,3)])
# plot_ts(RegFF, PP, bin_ts1, bin_ts2, "RegFF", "PP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test2 = bin_cor(RegBB, TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB")
binnedts      <- test2$Binned_time_series
bin_ts1 <- na.omit(test2$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test2$Binned_time_series[,c(1,3)])
plot_ts(RegBB, TEMP, bin_ts1, bin_ts2, "RegBB", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test2 = bin_cor(RegBB, TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB")
binnedts      <- test2$Binned_time_series
bin_ts1 <- na.omit(test2$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test2$Binned_time_series[,c(1,3)])
plot_ts(RegBB, TEMP, bin_ts1, bin_ts2, "RegBB", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test2 = bin_cor(RegBB, TEMP8000, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB")
binnedts      <- test2$Binned_time_series
bin_ts1 <- na.omit(test2$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test2$Binned_time_series[,c(1,3)])
plot_ts(RegBB, TEMP8000, bin_ts1, bin_ts2, "RegBB", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test2 = bin_cor(RegBB, PP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB")
# binnedts      <- test2$Binned_time_series
# bin_ts1 <- na.omit(test2$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test2$Binned_time_series[,c(1,3)])
# plot_ts(RegBB, PP, bin_ts1, bin_ts2, "RegBB", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test3 = bin_cor(FS, TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS")
binnedts      <- test3$Binned_time_series
bin_ts1 <- na.omit(test3$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test3$Binned_time_series[,c(1,3)])
plot_ts(FS, TEMP, bin_ts1, bin_ts2, "FS", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test3 = bin_cor(FS, TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS")
binnedts      <- test3$Binned_time_series
bin_ts1 <- na.omit(test3$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test3$Binned_time_series[,c(1,3)])
plot_ts(FS, TEMP, bin_ts1, bin_ts2, "FS", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test3 = bin_cor(FS, TEMP8000, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS")
binnedts      <- test3$Binned_time_series
bin_ts1 <- na.omit(test3$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test3$Binned_time_series[,c(1,3)])
plot_ts(FS, TEMP8000, bin_ts1, bin_ts2, "FS", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")
# 
# test3 = bin_cor(FS[1:8068,], TEMP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS")
# binnedts      <- test3$Binned_time_series
# bin_ts1 <- na.omit(test3$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test3$Binned_time_series[,c(1,3)])
# plot_ts(FS[1:8068,], TEMP, bin_ts1, bin_ts2, "FS", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")
# 
# test1 = bin_cor(FS[1:8068,], TEMP[1:8068,], FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF")
# binnedts<- test1$Binned_time_series
# bin_ts1 <- na.omit(test1$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test1$Binned_time_series[,c(1,3)])
# plot_ts(FS[1:8068,], TEMP[1:8068,], bin_ts1, bin_ts2, "RegFF", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")
# 
# 
# test3 = bin_cor(FS, PP, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_pp_FS")
# binnedts      <- test3$Binned_time_series
# bin_ts1 <- na.omit(test3$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test3$Binned_time_series[,c(1,3)])
# plot_ts(FS, PP, bin_ts1, bin_ts2, "FS", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

###Pollen et fire metrics

test4 = bin_cor(RegFF, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_PAR")
binnedts      <- test4$Binned_time_series
bin_ts1 <- na.omit(test4$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test4$Binned_time_series[,c(1,3)])
plot_ts(RegFF, PAR, bin_ts1, bin_ts2, "RegFF", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test4 = bin_cor(RegFF, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_PAR")
binnedts      <- test4$Binned_time_series
bin_ts1 <- na.omit(test4$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test4$Binned_time_series[,c(1,3)])
plot_ts(RegFF, PAR, bin_ts1, bin_ts2, "RegFF", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test4 = bin_cor(PP, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_PAR")
# binnedts      <- test4$Binned_time_series
# bin_ts1 <- na.omit(test4$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test4$Binned_time_series[,c(1,3)])
# plot_ts(PP, PAR, bin_ts1, bin_ts2, "RegFF", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test5 = bin_cor(RegBB, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_PAR")
binnedts      <- test5$Binned_time_series
bin_ts1 <- na.omit(test5$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test5$Binned_time_series[,c(1,3)])
plot_ts(RegBB, PAR, bin_ts1, bin_ts2, "RegBB", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test5 = bin_cor(RegBB, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_PAR")
binnedts      <- test5$Binned_time_series
bin_ts1 <- na.omit(test5$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test5$Binned_time_series[,c(1,3)])
plot_ts(RegBB, PAR, bin_ts1, bin_ts2, "RegBB", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test6 = bin_cor(FS, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
binnedts      <- test6$Binned_time_series
bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
plot_ts(FS, PAR, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test6 = bin_cor(FS, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
binnedts      <- test6$Binned_time_series
bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
plot_ts(FS, PAR, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], PAR, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")
# 
# test1 = bin_cor(PAR, TEMP[1:8068,], FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF")
# binnedts<- test1$Binned_time_series
# bin_ts1 <- na.omit(test1$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test1$Binned_time_series[,c(1,3)])
# plot_ts(PAR, TEMP[1:8068,], bin_ts1, bin_ts2, "RegFF", "TEMP", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

#Species

test7 = bin_cor(RegFF, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Picea")
binnedts      <- test7$Binned_time_series
bin_ts1 <- na.omit(test7$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test7$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Picea, bin_ts1, bin_ts2, "RegFF", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test7 = bin_cor(RegFF, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Picea")
binnedts      <- test7$Binned_time_series
bin_ts1 <- na.omit(test7$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test7$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Picea, bin_ts1, bin_ts2, "RegFF", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test8 = bin_cor(RegBB, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Picea")
binnedts      <- test8$Binned_time_series
bin_ts1 <- na.omit(test8$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test8$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Picea, bin_ts1, bin_ts2, "RegBB", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test8 = bin_cor(RegBB, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Picea")
binnedts      <- test8$Binned_time_series
bin_ts1 <- na.omit(test8$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test8$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Picea, bin_ts1, bin_ts2, "RegBB", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(FS, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(FS, Picea, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(FS, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(FS, Picea, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Picea, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Picea, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Picea, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Picea, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Picea, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test10 = bin_cor(RegFF, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Pinus")
binnedts      <- test10$Binned_time_series
bin_ts1 <- na.omit(test10$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test10$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pinus, bin_ts1, bin_ts2, "RegFF", "Pinus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test10 = bin_cor(RegFF, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Pinus")
binnedts      <- test10$Binned_time_series
bin_ts1 <- na.omit(test10$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test10$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pinus, bin_ts1, bin_ts2, "RegFF", "Pinus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test11 = bin_cor(RegBB, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Pinus")
binnedts      <- test11$Binned_time_series
bin_ts1 <- na.omit(test11$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test11$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pinus, bin_ts1, bin_ts2, "RegBB", "Pinus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test11 = bin_cor(RegBB, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Pinus")
binnedts      <- test11$Binned_time_series
bin_ts1 <- na.omit(test11$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test11$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pinus, bin_ts1, bin_ts2, "RegBB", "Pinus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test12 = bin_cor(FS, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Pinus")
binnedts      <- test12$Binned_time_series
bin_ts1 <- na.omit(test12$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test12$Binned_time_series[,c(1,3)])
plot_ts(FS, Pinus, bin_ts1, bin_ts2, "FS", "Pinus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test12 = bin_cor(FS, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Pinus")
binnedts      <- test12$Binned_time_series
bin_ts1 <- na.omit(test12$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test12$Binned_time_series[,c(1,3)])
plot_ts(FS, Pinus, bin_ts1, bin_ts2, "FS", "Pinus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Pinus, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Pinus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Pinus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Pinus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Pinus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test13 = bin_cor(RegFF, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test13$Binned_time_series
bin_ts1 <- na.omit(test13$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test13$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pinus, bin_ts1, bin_ts2, "RegFF", "Betula", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test13 = bin_cor(RegFF, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test13$Binned_time_series
bin_ts1 <- na.omit(test13$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test13$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pinus, bin_ts1, bin_ts2, "RegFF", "Betula", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test14 = bin_cor(RegBB, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020//2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test14$Binned_time_series
bin_ts1 <- na.omit(test14$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test14$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pinus, bin_ts1, bin_ts2, "RegBB", "Betula", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test14 = bin_cor(RegBB, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020//2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test14$Binned_time_series
bin_ts1 <- na.omit(test14$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test14$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pinus, bin_ts1, bin_ts2, "RegBB", "Betula", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test15 = bin_cor(FS, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test15$Binned_time_series
bin_ts1 <- na.omit(test15$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test15$Binned_time_series[,c(1,3)])
plot_ts(FS, Pinus, bin_ts1, bin_ts2, "FS", "Betula", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test15 = bin_cor(FS, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test15$Binned_time_series
bin_ts1 <- na.omit(test15$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test15$Binned_time_series[,c(1,3)])
plot_ts(FS, Pinus, bin_ts1, bin_ts2, "FS", "Betula", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Betula, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Betula, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Betula, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Betula, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Betula, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test16 = bin_cor(RegFF, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test16$Binned_time_series
bin_ts1 <- na.omit(test16$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test16$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pinus, bin_ts1, bin_ts2, "RegFF", "Alnus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test16 = bin_cor(RegFF, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test16$Binned_time_series
bin_ts1 <- na.omit(test16$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test16$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pinus, bin_ts1, bin_ts2, "RegFF", "Alnus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test17 = bin_cor(RegBB, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020//2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test17$Binned_time_series
bin_ts1 <- na.omit(test17$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test17$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pinus, bin_ts1, bin_ts2, "RegBB", "Alnus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test17 = bin_cor(RegBB, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020//2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test17$Binned_time_series
bin_ts1 <- na.omit(test17$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test17$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pinus, bin_ts1, bin_ts2, "RegBB", "Alnus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test18 = bin_cor(FS, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test18$Binned_time_series
bin_ts1 <- na.omit(test18$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test18$Binned_time_series[,c(1,3)])
plot_ts(FS, Pinus, bin_ts1, bin_ts2, "FS", "Alnus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test18 = bin_cor(FS, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test18$Binned_time_series
bin_ts1 <- na.omit(test18$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test18$Binned_time_series[,c(1,3)])
plot_ts(FS, Pinus, bin_ts1, bin_ts2, "FS", "Alnus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Alnus.c, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Alnus.c, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Alnus.c, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Alnus.c, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Alnus.c, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test19 = bin_cor(RegFF, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test19$Binned_time_series
bin_ts1 <- na.omit(test19$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test19$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Juniperus, bin_ts1, bin_ts2, "RegFF", "Juniperus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test19 = bin_cor(RegFF, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test19$Binned_time_series
bin_ts1 <- na.omit(test19$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test19$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Juniperus, bin_ts1, bin_ts2, "RegFF", "Juniperus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test20 = bin_cor(RegBB, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020//2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test20$Binned_time_series
bin_ts1 <- na.omit(test20$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test20$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Juniperus, bin_ts1, bin_ts2, "RegBB", "Juniperus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test20 = bin_cor(RegBB, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020//2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test20$Binned_time_series
bin_ts1 <- na.omit(test20$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test20$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Juniperus, bin_ts1, bin_ts2, "RegBB", "Juniperus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test21 = bin_cor(FS, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test21$Binned_time_series
bin_ts1 <- na.omit(test21$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test21$Binned_time_series[,c(1,3)])
plot_ts(FS, Juniperus, bin_ts1, bin_ts2, "FS", "Juniperus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test21 = bin_cor(FS, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test21$Binned_time_series
bin_ts1 <- na.omit(test21$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test21$Binned_time_series[,c(1,3)])
plot_ts(FS, Juniperus, bin_ts1, bin_ts2, "FS", "Juniperus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Juniperus, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Juniperus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Juniperus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Juniperus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Juniperus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test22 = bin_cor(RegFF, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test22$Binned_time_series
bin_ts1 <- na.omit(test22$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test22$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Populus, bin_ts1, bin_ts2, "RegFF", "Populus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test22 = bin_cor(RegFF, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test22$Binned_time_series
bin_ts1 <- na.omit(test22$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test22$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Populus, bin_ts1, bin_ts2, "RegFF", "Populus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test23 = bin_cor(RegBB, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test23$Binned_time_series
bin_ts1 <- na.omit(test23$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test23$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Populus, bin_ts1, bin_ts2, "RegBB", "Populus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test23 = bin_cor(RegBB, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test23$Binned_time_series
bin_ts1 <- na.omit(test23$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test23$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Populus, bin_ts1, bin_ts2, "RegBB", "Populus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test24 = bin_cor(FS, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test24$Binned_time_series
bin_ts1 <- na.omit(test24$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test24$Binned_time_series[,c(1,3)])
plot_ts(FS, Populus, bin_ts1, bin_ts2, "FS", "Populus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test24 = bin_cor(FS, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test24$Binned_time_series
bin_ts1 <- na.omit(test24$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test24$Binned_time_series[,c(1,3)])
plot_ts(FS, Populus, bin_ts1, bin_ts2, "FS", "Populus", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Populus, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Populus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Populus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Populus, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Populus, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test25 = bin_cor(RegFF, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test25$Binned_time_series
bin_ts1 <- na.omit(test25$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test25$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Myrica, bin_ts1, bin_ts2, "RegFF", "Myrica", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test25 = bin_cor(RegFF, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test25$Binned_time_series
bin_ts1 <- na.omit(test25$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test25$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Myrica, bin_ts1, bin_ts2, "RegFF", "Myrica", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test26 = bin_cor(RegBB, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test26$Binned_time_series
bin_ts1 <- na.omit(test26$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test26$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Myrica, bin_ts1, bin_ts2, "RegBB", "Myrica", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test26 = bin_cor(RegBB, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test26$Binned_time_series
bin_ts1 <- na.omit(test26$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test26$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Myrica, bin_ts1, bin_ts2, "RegBB", "Myrica", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test27 = bin_cor(FS, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test27$Binned_time_series
bin_ts1 <- na.omit(test27$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test27$Binned_time_series[,c(1,3)])
plot_ts(FS, Myrica, bin_ts1, bin_ts2, "FS", "Myrica", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test27 = bin_cor(FS, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test27$Binned_time_series
bin_ts1 <- na.omit(test27$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test27$Binned_time_series[,c(1,3)])
plot_ts(FS, Myrica, bin_ts1, bin_ts2, "FS", "Myrica", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Myrica, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Myrica, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Myrica, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Myrica, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Myrica, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test28 = bin_cor(RegFF, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test28$Binned_time_series
bin_ts1 <- na.omit(test28$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test28$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Lycopodium, bin_ts1, bin_ts2, "RegFF", "Lycopodium", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test28 = bin_cor(RegFF, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test28$Binned_time_series
bin_ts1 <- na.omit(test28$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test28$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Lycopodium, bin_ts1, bin_ts2, "RegFF", "Lycopodium", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test29 = bin_cor(RegBB, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test29$Binned_time_series
bin_ts1 <- na.omit(test29$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test29$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Lycopodium, bin_ts1, bin_ts2, "RegBB", "Lycopodium", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test29 = bin_cor(RegBB, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test29$Binned_time_series
bin_ts1 <- na.omit(test29$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test29$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Lycopodium, bin_ts1, bin_ts2, "RegBB", "Lycopodium", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test30 = bin_cor(FS, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test30$Binned_time_series
bin_ts1 <- na.omit(test30$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test30$Binned_time_series[,c(1,3)])
plot_ts(FS, Lycopodium, bin_ts1, bin_ts2, "FS", "Lycopodium", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test30 = bin_cor(FS, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test30$Binned_time_series
bin_ts1 <- na.omit(test30$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test30$Binned_time_series[,c(1,3)])
plot_ts(FS, Lycopodium, bin_ts1, bin_ts2, "FS", "Lycopodium", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Lycopodium, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Lycopodium, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Lycopodium, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Lycopodium, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Lycopodium, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test31 = bin_cor(RegFF, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test31$Binned_time_series
bin_ts1 <- na.omit(test31$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test31$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Salix, bin_ts1, bin_ts2, "RegFF", "Salix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test31 = bin_cor(RegFF, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test31$Binned_time_series
bin_ts1 <- na.omit(test31$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test31$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Salix, bin_ts1, bin_ts2, "RegFF", "Salix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test32 = bin_cor(RegBB, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test32$Binned_time_series
bin_ts1 <- na.omit(test32$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test32$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Salix, bin_ts1, bin_ts2, "RegBB", "Salix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test32 = bin_cor(RegBB, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test32$Binned_time_series
bin_ts1 <- na.omit(test32$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test32$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Salix, bin_ts1, bin_ts2, "RegBB", "Salix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test33 = bin_cor(FS, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test33$Binned_time_series
bin_ts1 <- na.omit(test33$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test33$Binned_time_series[,c(1,3)])
plot_ts(FS, Salix, bin_ts1, bin_ts2, "FS", "Salix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test33 = bin_cor(FS, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test33$Binned_time_series
bin_ts1 <- na.omit(test33$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test33$Binned_time_series[,c(1,3)])
plot_ts(FS, Salix, bin_ts1, bin_ts2, "FS", "Salix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Salix, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Salix, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Salix, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Salix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Salix, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test34 = bin_cor(RegFF, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test34$Binned_time_series
bin_ts1 <- na.omit(test34$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test34$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Larix, bin_ts1, bin_ts2, "RegFF", "Larix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test34 = bin_cor(RegFF, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test34$Binned_time_series
bin_ts1 <- na.omit(test34$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test34$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Larix, bin_ts1, bin_ts2, "RegFF", "Larix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test35 = bin_cor(RegBB, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test35$Binned_time_series
bin_ts1 <- na.omit(test35$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test35$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Larix, bin_ts1, bin_ts2, "RegBB", "Larix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test35 = bin_cor(RegBB, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test35$Binned_time_series
bin_ts1 <- na.omit(test35$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test35$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Larix, bin_ts1, bin_ts2, "RegBB", "Larix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test36 = bin_cor(FS, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test36$Binned_time_series
bin_ts1 <- na.omit(test36$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test36$Binned_time_series[,c(1,3)])
plot_ts(FS, Larix, bin_ts1, bin_ts2, "FS", "Larix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test36 = bin_cor(FS, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test36$Binned_time_series
bin_ts1 <- na.omit(test36$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test36$Binned_time_series[,c(1,3)])
plot_ts(FS, Larix, bin_ts1, bin_ts2, "FS", "Larix", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Larix, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Larix, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Larix, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Larix, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Larix, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test37 = bin_cor(RegFF, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test37$Binned_time_series
bin_ts1 <- na.omit(test37$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test37$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Cyperaceae, bin_ts1, bin_ts2, "RegFF", "Cyperaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test37 = bin_cor(RegFF, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test37$Binned_time_series
bin_ts1 <- na.omit(test37$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test37$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Cyperaceae, bin_ts1, bin_ts2, "RegFF", "Cyperaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test38 = bin_cor(RegBB, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test38$Binned_time_series
bin_ts1 <- na.omit(test38$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test38$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Cyperaceae, bin_ts1, bin_ts2, "RegBB", "Cyperaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test38 = bin_cor(RegBB, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test38$Binned_time_series
bin_ts1 <- na.omit(test38$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test38$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Cyperaceae, bin_ts1, bin_ts2, "RegBB", "Cyperaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test39 = bin_cor(FS, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test39$Binned_time_series
bin_ts1 <- na.omit(test39$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test39$Binned_time_series[,c(1,3)])
plot_ts(FS, Cyperaceae, bin_ts1, bin_ts2, "FS", "Cyperaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test39 = bin_cor(FS, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test39$Binned_time_series
bin_ts1 <- na.omit(test39$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test39$Binned_time_series[,c(1,3)])
plot_ts(FS, Cyperaceae, bin_ts1, bin_ts2, "FS", "Cyperaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Cyperaceae, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Cyperaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Cyperaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Cyperaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Cyperaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test40 = bin_cor(RegFF, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test40$Binned_time_series
bin_ts1 <- na.omit(test40$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test40$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Ericaceae, bin_ts1, bin_ts2, "RegFF", "Ericaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test40 = bin_cor(RegFF, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test40$Binned_time_series
bin_ts1 <- na.omit(test40$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test40$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Ericaceae, bin_ts1, bin_ts2, "RegFF", "Ericaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test41 = bin_cor(RegBB, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test41$Binned_time_series
bin_ts1 <- na.omit(test41$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test41$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Ericaceae, bin_ts1, bin_ts2, "RegBB", "Ericaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test41 = bin_cor(RegBB, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test41$Binned_time_series
bin_ts1 <- na.omit(test41$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test41$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Ericaceae, bin_ts1, bin_ts2, "RegBB", "Ericaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test42 = bin_cor(FS, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test42$Binned_time_series
bin_ts1 <- na.omit(test42$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test42$Binned_time_series[,c(1,3)])
plot_ts(FS, Ericaceae, bin_ts1, bin_ts2, "FS", "Ericaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test42 = bin_cor(FS, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test42$Binned_time_series
bin_ts1 <- na.omit(test42$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test42$Binned_time_series[,c(1,3)])
plot_ts(FS, Ericaceae, bin_ts1, bin_ts2, "FS", "Ericaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Ericaceae, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Ericaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Ericaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Ericaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Ericaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test43 = bin_cor(RegFF, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test43$Binned_time_series
bin_ts1 <- na.omit(test43$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test43$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Artemisia, bin_ts1, bin_ts2, "RegFF", "Artemisia", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test43 = bin_cor(RegFF, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test43$Binned_time_series
bin_ts1 <- na.omit(test43$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test43$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Artemisia, bin_ts1, bin_ts2, "RegFF", "Artemisia", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test44 = bin_cor(RegBB, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test44$Binned_time_series
bin_ts1 <- na.omit(test44$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test44$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Artemisia, bin_ts1, bin_ts2, "RegBB", "Artemisia", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test44 = bin_cor(RegBB, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test44$Binned_time_series
bin_ts1 <- na.omit(test44$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test44$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Artemisia, bin_ts1, bin_ts2, "RegBB", "Artemisia", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test45 = bin_cor(FS, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test45$Binned_time_series
bin_ts1 <- na.omit(test45$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test45$Binned_time_series[,c(1,3)])
plot_ts(FS, Artemisia, bin_ts1, bin_ts2, "FS", "Artemisia", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test45 = bin_cor(FS, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test45$Binned_time_series
bin_ts1 <- na.omit(test45$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test45$Binned_time_series[,c(1,3)])
plot_ts(FS, Artemisia, bin_ts1, bin_ts2, "FS", "Artemisia", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Artemisia, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Artemisia, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Artemisia, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Artemisia, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Artemisia, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test46 = bin_cor(RegFF, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test46$Binned_time_series
bin_ts1 <- na.omit(test46$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test46$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Poaceae, bin_ts1, bin_ts2, "RegFF", "Poaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test46 = bin_cor(RegFF, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test46$Binned_time_series
bin_ts1 <- na.omit(test46$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test46$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Poaceae, bin_ts1, bin_ts2, "RegFF", "Poaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test47 = bin_cor(RegBB, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test47$Binned_time_series
bin_ts1 <- na.omit(test47$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test47$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Poaceae, bin_ts1, bin_ts2, "RegBB", "Poaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test47 = bin_cor(RegBB, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test47$Binned_time_series
bin_ts1 <- na.omit(test47$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test47$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Poaceae, bin_ts1, bin_ts2, "RegBB", "Poaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test48 = bin_cor(FS, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test48$Binned_time_series
bin_ts1 <- na.omit(test48$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test48$Binned_time_series[,c(1,3)])
plot_ts(FS, Poaceae, bin_ts1, bin_ts2, "FS", "Poaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test48 = bin_cor(FS, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test48$Binned_time_series
bin_ts1 <- na.omit(test48$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test48$Binned_time_series[,c(1,3)])
plot_ts(FS, Poaceae, bin_ts1, bin_ts2, "FS", "Poaceae", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Poaceae, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Poaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Poaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Poaceae, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Poaceae, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test49 = bin_cor(RegFF, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test49$Binned_time_series
bin_ts1 <- na.omit(test49$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test49$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pediastrum, bin_ts1, bin_ts2, "RegFF", "Pediastrum", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test49 = bin_cor(RegFF, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test49$Binned_time_series
bin_ts1 <- na.omit(test49$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test49$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Pediastrum, bin_ts1, bin_ts2, "RegFF", "Pediastrum", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test50 = bin_cor(RegBB, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test50$Binned_time_series
bin_ts1 <- na.omit(test50$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test50$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pediastrum, bin_ts1, bin_ts2, "RegBB", "Pediastrum", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test50 = bin_cor(RegBB, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test50$Binned_time_series
bin_ts1 <- na.omit(test50$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test50$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Pediastrum, bin_ts1, bin_ts2, "RegBB", "Pediastrum", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test51 = bin_cor(FS, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test51$Binned_time_series
bin_ts1 <- na.omit(test51$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test51$Binned_time_series[,c(1,3)])
plot_ts(FS, Pediastrum, bin_ts1, bin_ts2, "FS", "Pediastrum", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test51 = bin_cor(FS, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test51$Binned_time_series
bin_ts1 <- na.omit(test51$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test51$Binned_time_series[,c(1,3)])
plot_ts(FS, Pediastrum, bin_ts1, bin_ts2, "FS", "Pediastrum", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Pediastrum, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Pediastrum, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Pediastrum, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Pediastrum, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Pediastrum, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test52 = bin_cor(RegFF, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test52$Binned_time_series
bin_ts1 <- na.omit(test52$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test52$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Potamogeton, bin_ts1, bin_ts2, "RegFF", "Potamogeton", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test52 = bin_cor(RegFF, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test52$Binned_time_series
bin_ts1 <- na.omit(test52$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test52$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Potamogeton, bin_ts1, bin_ts2, "RegFF", "Potamogeton", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test53 = bin_cor(RegBB, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test53$Binned_time_series
bin_ts1 <- na.omit(test53$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test53$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Potamogeton, bin_ts1, bin_ts2, "RegBB", "Potamogeton", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test53 = bin_cor(RegBB, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test53$Binned_time_series
bin_ts1 <- na.omit(test53$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test53$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Potamogeton, bin_ts1, bin_ts2, "RegBB", "Potamogeton", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test54 = bin_cor(FS, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test54$Binned_time_series
bin_ts1 <- na.omit(test54$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test54$Binned_time_series[,c(1,3)])
plot_ts(FS, Potamogeton, bin_ts1, bin_ts2, "FS", "Potamogeton", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test54 = bin_cor(FS, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test54$Binned_time_series
bin_ts1 <- na.omit(test54$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test54$Binned_time_series[,c(1,3)])
plot_ts(FS, Potamogeton, bin_ts1, bin_ts2, "FS", "Potamogeton", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Potamogeton, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Potamogeton, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Potamogeton, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Potamogeton, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Potamogeton, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test55 = bin_cor(RegFF, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test55$Binned_time_series
bin_ts1 <- na.omit(test55$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test55$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Nuphar, bin_ts1, bin_ts2, "RegFF", "Nuphar", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test55 = bin_cor(RegFF, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test55$Binned_time_series
bin_ts1 <- na.omit(test55$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test55$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Nuphar, bin_ts1, bin_ts2, "RegFF", "Nuphar", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test56 = bin_cor(RegBB, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test56$Binned_time_series
bin_ts1 <- na.omit(test56$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test56$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Nuphar, bin_ts1, bin_ts2, "RegBB", "Nuphar", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test56 = bin_cor(RegBB, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test56$Binned_time_series
bin_ts1 <- na.omit(test56$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test56$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Nuphar, bin_ts1, bin_ts2, "RegBB", "Nuphar", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test57 = bin_cor(FS, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test57$Binned_time_series
bin_ts1 <- na.omit(test57$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test57$Binned_time_series[,c(1,3)])
plot_ts(FS, Nuphar, bin_ts1, bin_ts2, "FS", "Nuphar", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test57 = bin_cor(FS, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test57$Binned_time_series
bin_ts1 <- na.omit(test57$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test57$Binned_time_series[,c(1,3)])
plot_ts(FS, Nuphar, bin_ts1, bin_ts2, "FS", "Nuphar", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Nuphar, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Nuphar, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Nuphar, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Nuphar, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Nuphar, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test58 = bin_cor(RegFF, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test58$Binned_time_series
bin_ts1 <- na.omit(test58$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test58$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Alnus.r, bin_ts1, bin_ts2, "RegFF", "Alnus.r", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test58 = bin_cor(RegFF, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegFF_Betula")
binnedts      <- test58$Binned_time_series
bin_ts1 <- na.omit(test58$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test58$Binned_time_series[,c(1,3)])
plot_ts(RegFF, Alnus.r, bin_ts1, bin_ts2, "RegFF", "Alnus.r", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test59 = bin_cor(RegBB, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test59$Binned_time_series
bin_ts1 <- na.omit(test59$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test59$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Alnus.r, bin_ts1, bin_ts2, "RegBB", "Alnus.r", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test59 = bin_cor(RegBB, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_RegBB_Betula")
binnedts      <- test59$Binned_time_series
bin_ts1 <- na.omit(test59$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test59$Binned_time_series[,c(1,3)])
plot_ts(RegBB, Alnus.r, bin_ts1, bin_ts2, "RegBB", "Alnus.r", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test60 = bin_cor(FS, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test60$Binned_time_series
bin_ts1 <- na.omit(test60$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test60$Binned_time_series[,c(1,3)])
plot_ts(FS, Alnus.r, bin_ts1, bin_ts2, "FS", "Alnus.r", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test60 = bin_cor(FS, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_Betula")
binnedts      <- test60$Binned_time_series
bin_ts1 <- na.omit(test60$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test60$Binned_time_series[,c(1,3)])
plot_ts(FS, Alnus.r, bin_ts1, bin_ts2, "FS", "Alnus.r", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# test6 = bin_cor(FS[1:8068,], Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_PAR")
# binnedts      <- test6$Binned_time_series
# bin_ts1 <- na.omit(test6$Binned_time_series[,1:2])
# bin_ts2 <- na.omit(test6$Binned_time_series[,c(1,3)])
# plot_ts(FS[1,8068,], Alnus.r, bin_ts1, bin_ts2, "FS", "PAR", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
# cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Alnus.r, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test9 = bin_cor(TEMP, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP, Alnus.r, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test9 = bin_cor(TEMP8000, Alnus.r, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_FS_picea")
binnedts      <- test9$Binned_time_series
bin_ts1 <- na.omit(test9$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test9$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, Alnus.r, bin_ts1, bin_ts2, "FS", "Picea", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

#TEMP vs PAR

test61 = bin_cor(TEMP, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_PAR")
binnedts<- test61$Binned_time_series
bin_ts1 <- na.omit(test61$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test61$Binned_time_series[,c(1,3)])
plot_ts(TEMP, PAR, bin_ts1, bin_ts2, "TEMP", "PAR", colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="y", KoCM="pearson")

test61 = bin_cor(TEMP, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_PAR")
binnedts<- test61$Binned_time_series
bin_ts1 <- na.omit(test61$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test61$Binned_time_series[,c(1,3)])
plot_ts(TEMP, PAR, bin_ts1, bin_ts2, "TEMP", "PAR", colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

test61 = bin_cor(TEMP8000, PAR, FLAGTAU=3, ofilename="E:/Sauvegarde_24052020/2-En_cours/AXE_2_FEUX_PASSES/Bincor/ccf_temp_PAR")
binnedts<- test61$Binned_time_series
bin_ts1 <- na.omit(test61$Binned_time_series[,1:2])
bin_ts2 <- na.omit(test61$Binned_time_series[,c(1,3)])
plot_ts(TEMP8000, PAR, bin_ts1, bin_ts2, "TEMP", "PAR", colts2=2, colbints1=3, colbints2=4, device="screen")
cor_ts(bin_ts1, bin_ts2, rmltrd="n", KoCM="pearson")

# Sedimentation rate

library(ggplot2)

ggplot(data = sedim, aes(x = ï..Age, y = Taux.sÃ.dimentation)) +
  geom_area(fill = "black") +
  theme_classic()
