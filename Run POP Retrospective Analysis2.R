###################################################################################
# Script for finalizing GOA rockfish models
# First section of code adapted from P. Hulson retrospective code used in 2011
# 2014 revision by D. Hanselman
# Generalized code to work directly from the .dat file for each rockfish
# Requires a models.dat file to read .ctls in for sensitivity, and for base model ctl file
# Runs retrospective models, will also run MCMC (set flag mcmcon<-"YES")
# Makes copies of results files (and MCMC evalout when MCMC is on)
# Runs sensitivity models for sensitivity graph
# Plots retrospective plots
# Plots some posterior and prior distribution plots for M and q
###################################################################################

################ Load up some libraries needed
library(ggplot2)
library(lattice)
library(reshape2)
library(MASS)
library(emdbook)


#============================================================================================
#============= Set up directories
#============================================================================================
setwd("C:/SA/Retro")
path<-getwd()

CTL<-read.table(paste(path,"/POPNov2011/models.dat",sep=""),sep="\t")

path
pathD<-paste(path,"/Data",sep="")

pathR<-paste(path,"/Results",sep="")

pathM<-paste(path,"/Model",sep="")

# Get current data file

DAT<-readLines(paste(pathD,"/goa_pop_2011.dat",sep=""),warn=FALSE)

Sec_st<-grep("#-",DAT)
Sec_end<-grep("#!",DAT)

st_end<-matrix(NA,nrow=length(Sec_st),ncol=2)
st_end[,1]<-Sec_st
st_end[,2]<-Sec_end

mcmcon<-"NO"
mcmcruns<-50000  # Could change these, but I like 5000 as a manageable number to deal with
mcmcsave<-mcmcruns/5000
#============================================================================================
#============= Set up some model dimensions
#============================================================================================
styr<-as.numeric(DAT[Sec_st[2]-3]) # start of model (example 1961 for POP)
modelyear<-as.numeric(DAT[Sec_st[2]-1]) #current model year
nages<-as.numeric(DAT[Sec_st[2]+3]) # number of age bins
nlens<-as.numeric(DAT[Sec_st[2]+5]) # number of length bins
numretros<-10 # number of retrospective years

#============================================================================================
#============= Run retrospective loop
#============================================================================================
T_start<-Sys.time() Timer start

for(y in 1:(numretros+1)) {

#============================================================================================
#============= Concatenate DAT file
#============================================================================================

# Set endyr
yrs_retro<-seq(modelyear-numretros,modelyear)
endyr<-yrs_retro[numretros-y+2]
nyrs<-endyr-styr+1
DAT_retro<-c(DAT[st_end[1,1]:st_end[1,2]],as.character(endyr),DAT[st_end[2,1]:st_end[2,2]])
# This is a new branch
# Fishery catch
DAT_retro<-c(DAT_retro,paste(scan(text=DAT[Sec_st[3]-1])[1:nyrs],collapse=" "),DAT[st_end[3,1]:st_end[3,2]])
# Trawl survey biomass
 BTSb_yrs<-length(which(scan(text=DAT[Sec_st[5]-1])<(endyr+1)))
DAT_retro<-c(
DAT_retro,
as.character(BTSb_yrs),
DAT[st_end[4,1]:st_end[4,2]],
paste(scan(text=DAT[Sec_st[5]-1])[1:BTSb_yrs],collapse=" "),
DAT[st_end[5,1]:st_end[5,2]],
paste(scan(text=DAT[Sec_st[6]-1])[1:BTSb_yrs],collapse=" "),
DAT[st_end[6,1]:st_end[6,2]],
paste(scan(text=DAT[Sec_st[7]-1])[1:BTSb_yrs],collapse=" "),
DAT[st_end[7,1]:st_end[7,2]],
paste(scan(text=DAT[Sec_st[8]-1])[1:BTSb_yrs],collapse=" "),
DAT[st_end[8,1]:st_end[8,2]],
paste(scan(text=DAT[Sec_st[9]-1])[1:BTSb_yrs],collapse=" "),
DAT[st_end[9,1]:st_end[9,2]])

# Fish age comp
FAC_yrs<-length(which(scan(text=DAT[Sec_st[11]-1])<(endyr)))
DAT_retro<-c(DAT_retro,
as.character(FAC_yrs),
DAT[st_end[10,1]:st_end[10,2]],
paste(scan(text=DAT[Sec_st[11]-1])[1:FAC_yrs],collapse=" "),
DAT[st_end[11,1]:st_end[11,2]],
paste(scan(text=DAT[Sec_st[12]-1])[1:FAC_yrs],collapse=" "),
DAT[st_end[12,1]:st_end[12,2]],
paste(scan(text=DAT[Sec_st[13]-1])[1:FAC_yrs],collapse=" "),
DAT[st_end[13,1]:st_end[13,2]],
paste(scan(text=DAT[Sec_st[14]-1])[1:FAC_yrs],collapse=" "),
DAT[st_end[14,1]:st_end[14,2]])
for(i in 1:FAC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=DAT[Sec_st[15]-FAC_yrs-1+i]) ,collapse = " "))

DAT_retro<-c(DAT_retro,DAT[st_end[15,1]:st_end[15,2]])

# Survey age comp
 SAC_yrs<-length(which(scan(text=DAT[Sec_st[17]-1])<(endyr)))
 DAT_retro<-c(DAT_retro,
  as.character(SAC_yrs),
  DAT[st_end[16,1]:st_end[16,2]],
  paste(scan(text=DAT[Sec_st[17]-1])[1:SAC_yrs],collapse=" "),
  DAT[st_end[17,1]:st_end[17,2]], DAT[st_end[11,1]:st_end[11,2]],
  paste(scan(text=DAT[Sec_st[18]-1])[1:SAC_yrs],collapse=" "),
  DAT[st_end[18,1]:st_end[18,2]],
  paste(scan(text=DAT[Sec_st[19]-1])[1:SAC_yrs],collapse=" "),
  DAT[st_end[19,1]:st_end[19,2]],
  paste(scan(text=DAT[Sec_st[20]-1])[1:SAC_yrs],collapse=" "),
  DAT[st_end[20,1]:st_end[20,2]])
  for(i in 1:SAC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=DAT[Sec_st[21]-SAC_yrs-1+i]) ,collapse = " "))

 DAT_retro<-c(DAT_retro,DAT[st_end[21,1]:st_end[21,2]])


# Fish size comp
 FSC_yrs<-length(which(scan(text=DAT[Sec_st[23]-1])<(endyr)))
 
 DAT_retro<-c(DAT_retro,
  as.character(FSC_yrs),
  DAT[st_end[22,1]:st_end[22,2]],
  paste(scan(text=DAT[Sec_st[23]-1])[1:FSC_yrs],collapse=" "),
  DAT[st_end[23,1]:st_end[23,2]],
  paste(scan(text=DAT[Sec_st[24]-1])[1:FSC_yrs],collapse=" "),
  DAT[st_end[24,1]:st_end[24,2]],
  paste(scan(text=DAT[Sec_st[25]-1])[1:FSC_yrs],collapse=" "),
  DAT[st_end[25,1]:st_end[25,2]],
  paste(scan(text=DAT[Sec_st[26]-1])[1:FSC_yrs],collapse=" "),
  DAT[st_end[26,1]:st_end[26,2]])
  for(i in 1:FSC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=DAT[Sec_st[27]-FSC_yrs-1+i]) ,collapse = " "))

DAT_retro<-c(DAT_retro,DAT[st_end[27,1]:st_end[27,2]])

# Survey size comp
 SSC_yrs<-length(which(scan(text=DAT[Sec_st[29]-1])<(endyr+1)))

DAT_retro<-c(DAT_retro,
 as.character(SSC_yrs),
 DAT[st_end[28,1]:st_end[28,2]],
 paste(scan(text=DAT[Sec_st[29]-1])[1:SSC_yrs],collapse=" "),
 DAT[st_end[29,1]:st_end[29,2]],
 paste(scan(text=DAT[Sec_st[30]-1])[1:SSC_yrs],collapse=" "),
 DAT[st_end[30,1]:st_end[30,2]],
 paste(scan(text=DAT[Sec_st[31]-1])[1:SSC_yrs],collapse=" "),
 DAT[st_end[31,1]:st_end[31,2]],
 paste(scan(text=DAT[Sec_st[32]-1])[1:SSC_yrs],collapse=" "),
 DAT[st_end[32,1]:st_end[32,2]])
 for(i in 1:SSC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=DAT[Sec_st[33]-SSC_yrs-1+i]) ,collapse = " "))

DAT_retro<-c(DAT_retro,DAT[st_end[33,1]:st_end[33,2]])

# Write data file
write.table(DAT_retro,file=paste(pathM,"/goa_pop_2011.dat",sep=""),quote=F,row.names=F,col.names=F)


#============================================================================================
#============= Run model
#============================================================================================

## set your number of MCMC runs at the top of the program... 
setwd(pathM)
if(mcmcon=="YES") shell(paste('tem.EXE',' -mcmc ',mcmcruns,'-mcsave ',mcmcsave)) else
   shell(paste('tem.exe ','-nox'))
   


#============================================================================================
#============= Get/write results
#============================================================================================


if(mcmcon=="YES") {
  shell(paste('"tem.exe"', '-mceval'))
  file.copy(from=paste(pathM,"/evalout.prj",sep=""),to=paste(pathR,"/mcmc_",modelyear-(y-1),".std",sep=""),overwrite=T)
 }

file.copy(from=paste(pathM,"/tem.STD",sep=""),to=paste(pathR,"/std_",modelyear-(y-1),".std",sep=""),overwrite=T)
file.copy(from=paste(pathM,"/report.rep",sep=""),to=paste(pathR,"/rep_",modelyear-(y-1),".rep",sep=""),overwrite=T)
file.copy(from=paste(pathM,"/tem.par",sep=""),to=paste(pathR,"/par_",modelyear-(y-1),".par",sep=""),overwrite=T)
file.copy(from=paste(pathM,"/proj.dat",sep=""),to=paste(pathR,"/prj_",modelyear-(y-1),".prj",sep=""),overwrite=T)


}

T_end<-Sys.time()

T_end-T_start
#---------------------------------------------
# End of retrospective model running loop
#---------------------------------------------
#
#-------------------------
#Reset data file to run sensitivities
#-------------------------

setwd(pathM)
#write.table(DAT,file=paste(pathM,"/goa_pop_2011.dat",sep=""),,quote=F,,,,,row.names=F,col.names=F)
write.table(DAT,file=paste(pathM,"/goa_pop_2011.dat",sep=""),quote=F,row.names=F,col.names=F)

#############################
# Run sensitivity runs
##############################
for(y in 1:length(CTL[1,])) {
  write.table(CTL[,y],file=paste(pathM,"/tem.ctl",sep=""),quote=F,row.names=F,col.names=F)
  shell(paste('tem.exe ','-nox'))

#============================================================================================
#============= Get/write results
#============================================================================================


file.copy(from=paste(pathM,"/tem.STD",sep=""),to=paste(pathR,"/std_sens_",y,".std",sep=""),overwrite=T)
file.copy(from=paste(pathM,"/report.rep",sep=""),to=paste(pathR,"/rep_sens_",y,".rep",sep=""),overwrite=T)
file.copy(from=paste(pathM,"/tem.par",sep=""),to=paste(pathR,"/par_sens_",y,".par",sep=""),overwrite=T)
file.copy(from=paste(pathM,"/proj.dat",sep=""),to=paste(pathR,"/par_prj_",y,".prj",sep=""),overwrite=T)

}
#============================================================================================
#============= Plot sensitivity graph
#============================================================================================


### Little function for extracting stuff out of the sensitivity runs..
filesearch <- function (x,y) {
  output<-cbind(seq(1,y),seq(1,y),seq(1,y),seq(1,y),seq(1,y),seq(1,y),seq(1,y))
  colnames(output)<-c("Test","lnL","ABC","SSB","Catchability","Natmort","MeanRec")
  for (i in 1:y) {
    z<-paste(x,i,".rep",sep="")
    f <- readLines(z)
    modname<-f[grep("Model name",f,
                    value=FALSE)+1]
    ABC<- f[grep(" ABC for 2012",f,
                 value=FALSE)+1]
    Like<-f[grep("Total likelihood",f,
                 value=FALSE)+1]
    SSB<-f[grep("Female_Spawning Biomass for 2012",f,
                value=FALSE)+1] 
    q<-f[grep("q_trawl",f,value=FALSE)+1]
    M<-f[grep("nat_mort",f,value=FALSE)+1]
    rec<-f[grep("log_mean_rec",f,value=FALSE)+1]
    
    output[i,1]<-modname
    output[i,2] <-as.numeric(Like)
    output[i,3] <- as.numeric(ABC)
    output[i,4] <- as.numeric(SSB)
    output[i,5] <- as.numeric(strsplit(q," ")[[1]])[1]
    output[i,6] <- as.numeric(strsplit(M," ")[[1]])[1]
    output[i,7] <- as.numeric(strsplit(rec," ")[[1]])[1]
    
  }
  return(output)
}

#location of data files
setwd(pathR)
y<-filesearch("rep_sens_",21)
#sensplot<-function(x) {
y<-data.frame(y,stringsAsFactors = FALSE)
for(j in 3:7) {
  y[,j]<-(as.numeric(y[,j])/as.numeric(y[1,j])*100)-100 }
y[,2]<-(as.numeric(y[,2])-as.numeric(y[1,2]))
z<-melt(y[-1,],id=c("Test"))
z[["sign"]]=ifelse(z[["value"]]>=0,"positive","negative")
a<-  ggplot(data=z, aes(x=Test, y=value,fill=sign))+scale_fill_manual(values = c("positive" = "darkblue", "negative" = "red")) 
b <- a + geom_bar(stat = "identity", position = "stack")
b <- b + scale_fill_brewer(palette = "Set1")
sensitivity_theme <- theme_update(axis.text.x = element_text(angle = 90,
                                                             hjust = 1), panel.grid.major = element_line(colour = "grey90"),
                                  panel.grid.minor = element_blank(), panel.background = element_blank(),
                                  axis.ticks = element_blank(), legend.position = "none")
c <- b + facet_grid(variable ~ .) + theme(legend.position = "none")
c2 <- c + facet_grid(variable ~ ., scale = "free_y")
c2  <-c2+labs(y="Percent difference from reference model")+geom_abline(intercept = 0, ,slope=0,colour = "dark green", size = 1.)
c2

#============================================================================================
#============= Plot retro graphs
#============================================================================================

  setwd(pathR)
  ssb<-seq(1,(modelyear-styr+1))
  totbio<-seq(1,(modelyear-styr+1))
  i<-0
  for (i in modelyear:(modelyear-numretros)) {
   print(i)
    z<-paste("rep_",i,".rep",sep="")
    f <- readLines(z)
    ssb1<-f[grep("SpBiom",f,value=FALSE)]
    ssb1<-sub("SpBiom  ","",ssb1)
    ssb1<-scan(text=ssb1)
    length(ssb1)<-modelyear-styr+1
    ssb<-cbind(ssb,ssb1)
   tb1<-f[grep("Tot_biom",f,value=FALSE)]
   tb1<-sub("Tot_biom  ","",tb1)
   tb1<-scan(text=tb1)
   length(tb1)<-modelyear-styr+1
   totbio<-cbind(totbio,tb1)
  print(i)
  }
ssb[,1]<-seq((modelyear-length(ssb[,1])+1),modelyear)
colnames(ssb)<-c("Year",seq(modelyear,modelyear-numretros))
totbio[,1]<-seq((modelyear-length(totbio[,1])+1),modelyear)
colnames(totbio)<-c("Year",seq(modelyear,modelyear-numretros))

#### Do recruitment from std files

### whitespace trimmer helper function
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
###

recdevs<-seq(1,(modelyear-styr+1+nages-2))
i<-0
for (i in modelyear:(modelyear-numretros)) {
  print(i)
  z<-paste("std_",i,".std",sep="")
  f <- readLines(z)
  g<-trim(f)
  xx<-strsplit(g[grep("log_rec_dev",g,value=FALSE)]," ")
  xx<-as.numeric(sapply(xx, "[", 9))
  length(xx)<-modelyear-styr+1+nages-2
  recdevs<-cbind(recdevs,xx)
    print(i)
}

 recdevs[,1]<-seq((modelyear-length(recdevs[,1])+1
                   ),modelyear)
 colnames(recdevs)<-c("Year",seq(modelyear,modelyear-numretros))
                        
#########################################
# Do some retro plots
#########################################
## Simple brute force retro line graph

plot(ssb[,1],100*(ssb[,2]-ssb[,2])/ssb[,2],type="l",lwd=3,xaxt="n",las=2,xlab="",ylab="",cex.axis=1.5,ylim=c(-50,50),lty=2)

lines(ssb[,1],100*(ssb[,3]-ssb[,2])/ssb[,2],lwd=1.5,col="purple")
lines(ssb[,1],100*(ssb[,4]-ssb[,2])/ssb[,2],lwd=1.5,col="violet")
lines(ssb[,1],100*(ssb[,5]-ssb[,2])/ssb[,2],lwd=1.5,col="dark blue")
lines(ssb[,1],100*(ssb[,6]-ssb[,2])/ssb[,2],lwd=1.5,col="blue")
lines(ssb[,1],100*(ssb[,7]-ssb[,2])/ssb[,2],lwd=1.5,col="green")
lines(ssb[,1],100*(ssb[,8]-ssb[,2])/ssb[,2],lwd=1.5,col="dark green")
lines(ssb[,1],100*(ssb[,9]-ssb[,2])/ssb[,2],lwd=1.5,col="red")
lines(ssb[,1],100*(ssb[,10]-ssb[,2])/ssb[,2],lwd=1.5,col="dark red")
lines(ssb[,1],100*(ssb[,11]-ssb[,2])/ssb[,2],lwd=1.5,col="brown")

mtext("Percent differences",side=2,line=6,cex=1.1)
mtext("from terminal year",side=2,line=4.5,cex=1.1)

for(y in seq(1,length(ssb[,1]),2)){
  axis(side=1,at=ssb[,1][y],ssb[,1][y],cex.axis=1.5,las=2)}

mtext("Year",side=1,line=4.25,cex=1.1)

################ Total bio


plot(totbio[,1],100*(totbio[,2]-totbio[,2])/totbio[,2],type="l",lwd=3,xaxt="n",las=2,xlab="",ylab="",cex.axis=1.5,ylim=c(-50,50),lty=2)

lines(totbio[,1],100*(totbio[,3]-totbio[,2])/totbio[,2],lwd=1.5,col="purple")
lines(totbio[,1],100*(totbio[,4]-totbio[,2])/totbio[,2],lwd=1.5,col="violet")
lines(totbio[,1],100*(totbio[,5]-totbio[,2])/totbio[,2],lwd=1.5,col="dark blue")
lines(totbio[,1],100*(totbio[,6]-totbio[,2])/totbio[,2],lwd=1.5,col="blue")
lines(totbio[,1],100*(totbio[,7]-totbio[,2])/totbio[,2],lwd=1.5,col="green")
lines(totbio[,1],100*(totbio[,8]-totbio[,2])/totbio[,2],lwd=1.5,col="dark green")
lines(totbio[,1],100*(totbio[,9]-totbio[,2])/totbio[,2],lwd=1.5,col="red")
lines(totbio[,1],100*(totbio[,10]-totbio[,2])/totbio[,2],lwd=1.5,col="dark red")
lines(totbio[,1],100*(totbio[,11]-totbio[,2])/totbio[,2],lwd=1.5,col="brown")

mtext("Percent differences",side=2,line=6,cex=1.1)
mtext("from terminal year",side=2,line=4.5,cex=1.1)

for(y in seq(1,length(totbio[,1]),2)){
  axis(side=1,at=totbio[,1][y],totbio[,1][y],cex.axis=1.5,las=2)}

mtext("Year",side=1,line=4.25,cex=1.1)
################ Do cool retro plots
recbak<-recdevs
recdevs<-recdevs[-seq(1,41),]
plot(recdevs[,1],100*(recdevs[,2]-recdevs[,2])/recdevs[,2],type="l",lwd=3,xaxt="n",las=2,xlab="",ylab="",cex.axis=1.5,ylim=c(-150,150),lty=2)

lines(recdevs[,1],100*(recdevs[,3]-recdevs[,2])/recdevs[,2],lwd=1.5,col="purple")
lines(recdevs[,1],100*(recdevs[,4]-recdevs[,2])/recdevs[,2],lwd=1.5,col="violet")
lines(recdevs[,1],100*(recdevs[,5]-recdevs[,2])/recdevs[,2],lwd=1.5,col="dark blue")
lines(recdevs[,1],100*(recdevs[,6]-recdevs[,2])/recdevs[,2],lwd=1.5,col="blue")
lines(recdevs[,1],100*(recdevs[,7]-recdevs[,2])/recdevs[,2],lwd=1.5,col="green")
lines(recdevs[,1],100*(recdevs[,8]-recdevs[,2])/recdevs[,2],lwd=1.5,col="dark green")
lines(recdevs[,1],100*(recdevs[,9]-recdevs[,2])/recdevs[,2],lwd=1.5,col="red")
lines(recdevs[,1],100*(recdevs[,10]-recdevs[,2])/recdevs[,2],lwd=1.5,col="dark red")
lines(recdevs[,1],100*(recdevs[,11]-recdevs[,2])/recdevs[,2],lwd=1.5,col="brown")

mtext("Percent differences",side=2,line=6,cex=1.1)
mtext("from terminal year",side=2,line=4.5,cex=1.1)

for(y in seq(1,length(recdevs[,1]),2)){
  axis(side=1,at=recdevs[,1][y],recdevs[,1][y],cex.axis=1.5,las=2)}

mtext("Year",side=1,line=4.25,cex=1.1)

##### Begin squid plot for recruitment cohort changes over time
### helper function for adding shaded polygons
addpoly <- function(yrvec, lower, upper, shadecol = rgb(0, 
                                                        0, 0, 0.1), col = 1) {
  polygon(x = c(yrvec, rev(yrvec)), y = c(lower, rev(upper)), 
          border = NA, col = shadecol)
  lines(yrvec, lower, lty = 3, col = col)
  lines(yrvec, upper, lty = 3, col = col)
}
####
### Squids
#### Make some color palletes
colvec <- rainbow(numretros+1, alpha = 0.7)
shadecolvec <- rainbow(numretros+1, alpha = 0.075)
colvec.txt <- colvec
for (i in 1:length(colvec)) {
  tmp <- col2rgb(colvec[i])/255
  colvec.txt[i] <- rgb(tmp[1]/2, tmp[2]/2, tmp[3]/2, alpha = 0.7)
}
print(cbind(colvec, colvec.txt))
ylim <-c(-1., 1)
ylim <- ylim + c(0, 1)
xlim <- c(0, 10)
  xlim <- xlim + c(-0.8, 0.8)
labies<-seq(1997,2008)
plot(0, type = "n", xlim = xlim, ylim = ylim, ylab = "Recruitment deviation", 
     xlab = "Years since birth", axes = FALSE,main="Pacific ocean perch recruitment retrospective")
axis(1, at = 0:10)
axis(2, at = -1:1., las = 1)
abline(h = 0, col = "grey")
box()

  i=0
  for(iy in 32:22) {
    i=i+1
    lines(seq(0,10),c(rev(recdevs[iy,2:12])[(12-i):11],rep(NA,11-i)), 
          type = "o", col = colvec[iy-21], lwd = 3, pch = 16)
       text(x=i-1+0.5,y=recdevs[iy,2],labels=as.character(2012-i),col=colvec[iy-21]) }

  legend("top", lwd = 3, lty = 1, pch = 16, col = colvec, 
         legend = seq(2001,2011), title = "Cohort birth year", ncol = 6, 
         bg = rgb(1, 1, 1, 0.3), box.col = NA,cex=1)

### Next step, generalize the squid, move on to MCMC retrospective plots... sweet!

#### Read in mcmc files for retro:
#Example of whats in a evalout for 2011
#totbio and ssb  102
#catch and ssb projections	30
#params	9
#recdevs	73
#rec_proj	10
#next year tot bio	1
#total	225

(for i in modelyear:(modelyear-numretros))
mcmc<-(read.table("mcmc_2011.std",header=F,sep="")) #,nrow=5000,ncol=50+3*(i-styr+1)+nages-recage)
mcmc<-mcmc[1001:5000,] # remove burning

### Identify variables
mcmcnames<-  c("sigr","q1",  "q2",	"f40",	"M",	"ssb_next",	"ABC",	"obj_fun")
for(i in styr:endyr) mcmcnames<- c(mcmcnames,paste("totbio",i,sep="")) 
for(i in (styr-nages+2):endyr) mcmcnames<- c(mcmcnames,paste("recdev",i,sep="")) 
for(i in styr:endyr) mcmcnames<- c(mcmcnames,paste("ssb",i,sep="")) 
mcmcnames<-c(mcmcnames,"LMR")
for(i in (endyr+1):(endyr+15)) mcmcnames<-c(mcmcnames,paste("ssbproj",i,sep=""))
for(i in (endyr+1):(endyr+15)) mcmcnames<-c(mcmcnames,paste("catchproj",i,sep=""))
for(i in (endyr+1):(endyr+10)) mcmcnames<-c(mcmcnames,paste("recproj",i,sep=""))
mcmcnames<-c(mcmcnames,paste("totbioproj",endyr+1,sep=""))
names(mcmc)<-mcmcnames

### Do percentile confidence bounds
### setup vectors

mcmcmed<-mcmc[1,]
mcmcuci<-mcmc[1,]
mcmclci<-mcmc[1,]
mcmed<-mcmc[1,]
mcuci<-mcmc[1,]
mclci<-mcmc[1,]
ssbmed<-data.frame(matrix(NA,nrow=11,ncol=51))
ssblci<-data.frame(matrix(NA,nrow=11,ncol=51))
ssbuci<-data.frame(matrix(NA,nrow=11,ncol=51))
l<-0
for (j in modelyear:(modelyear-numretros)) {
l=l+1
  mcmc<-(read.table(paste("mcmc_",j,".std",sep=""),header=F,sep="")) #,nrow=5000,ncol=50+3*(i-styr+1)+nages-recage)
  mcmc<-mcmc[1001:5000,] # remove burning
  
  for (i in 1:length(mcmc[1,])) {
    
    mcmcmed[i]<-median(mcmc[,i])
    mcmclci[i]<-quantile(mcmc[,i],0.025)
    mcmcuci[i]<-quantile(mcmc[,i],0.975) }
 # mcmed<-rbind(mcmed,mcmcmed)
#  mcuci<-rbind(mcuci,mcmcuci)
#  mclci<-rbind(mclci,mcmclci)
  ssbmed[l,1:(endyr-styr+1-l+1)]<-mcmcmed[(133-2*(l-1)):(183-2*(l-1)-(l-1))]
  ssblci[l,1:(endyr-styr+1-l+1)]<-mcmclci[(133-2*(l-1)):(183-2*(l-1)-(l-1))]
  ssbuci[l,1:(endyr-styr+1-l+1)]<-mcmcuci[(133-2*(l-1)):(183-2*(l-1)-(l-1))]
}

plot(styr:modelyear,ssbuci[1,]/1000,pch="",ylim=c(0,(1.1*max(ssbuci,na.rm=T)/1000)),xlab="Year",ylab="Spawning biomass (kt)")
for (i in 1:numretros) {
  lines(styr:(modelyear-i+1),ssbmed[i,1:(modelyear-styr-i+2)]/1000,col=colvec[i],lwd=2)
  addpoly(styr:(modelyear-i+1),ssblci[i,1:(modelyear-styr-i+2)]/1000,ssbuci[i,1:(modelyear-styr-i+2)]/1000,shadecol=shadecolvec[i])
  
}

plot(styr:modelyear,ssbuci[1,]/1000,pch="",ylim=c(-100,100),xlab="Year",ylab="Percent difference from terminal year")
abline(h=0,lwd=3,lty=2,col=colvec[1])
for (i in 1:numretros) {
  lines(styr:(modelyear-i+1),(1-ssbmed[i,1:(modelyear-styr-i+2)]/(ssbmed[1,1:(modelyear-styr-i+2)]))*100,col=colvec[i])
  addpoly(styr:(modelyear-i+1),(1-ssblci[i,1:(modelyear-styr-i+2)]/(ssbmed[1,1:(modelyear-styr-i+2)]))*100,
          (1- ssbuci[i,1:(modelyear-styr-i+2)]/(ssbmed[1,1:(modelyear-styr-i+2)]))*100,shadecol=shadecolvec[i])
  
}

## Do some posterior/prior plots

mcmc<-(read.table(paste("mcmc_",modelyear,".std",sep=""),header=F,sep="")) #,nrow=5000,ncol=50+3*(i-styr+1)+nages-recage)
mcmc<-mcmc[1001:5000,] # remove burning

### Identify variables
mcmcnames<-  c("sigr","q1",  "q2",  "f40",	"M",	"ssb_next",	"ABC",	"obj_fun")
for(i in styr:endyr) mcmcnames<- c(mcmcnames,paste("totbio",i,sep="")) 
for(i in (styr-nages+2):endyr) mcmcnames<- c(mcmcnames,paste("recdev",i,sep="")) 
for(i in styr:endyr) mcmcnames<- c(mcmcnames,paste("ssb",i,sep="")) 
mcmcnames<-c(mcmcnames,"LMR")
for(i in (endyr+1):(endyr+15)) mcmcnames<-c(mcmcnames,paste("ssbproj",i,sep=""))
for(i in (endyr+1):(endyr+15)) mcmcnames<-c(mcmcnames,paste("catchproj",i,sep=""))
for(i in (endyr+1):(endyr+10)) mcmcnames<-c(mcmcnames,paste("recproj",i,sep=""))
mcmcnames<-c(mcmcnames,paste("totbioproj",endyr+1,sep=""))
names(mcmc)<-mcmcnames
par(bg="transparent",col.axis="goldenrod",col="goldenrod",col.main="white",col.lab="goldenrod")
# Filled Density Plot
# POP catchability
basectl<-scan(text=as.character(CTL[3:59,1]))
qpr<-rlnorm(40000,log(as.numeric(basectl[26])),as.numeric(basectl[27])) ## read values from CTL file

d <- density(qpr)
e<- density(mcmc$q1)
plot(d, main="Prior and posterior of q",xlim=c(0,5),ylim=c(0,1.1),xlab="Catchability",yaxt="n")
polygon(d, col="red", border="blue")
polygon(e, col=rgb(0.2,0.7 , 0.8,0.5), border="red")

# POP natural mortality
mpr<-rlnorm(40000,log(as.numeric(basectl[17])),as.numeric(basectl[18])) ## read values from CTL file

d <- density(mpr)
e<- density(mcmc$M)
plot(d, main="Prior and posterior of q",xlim=c(0,0.1),ylim=c(0,100),xlab="Natural mortality",yaxt="n")
polygon(d, col="red", border="blue")
polygon(e, col=rgb(0.2,0.7 , 0.8,0.5), border="red")



z<-cbind(mcmc$M,mcmc$q1)
HPDregionplot(mcmc(z),col="green",lwd=3,xlab="Natural mortality",ylab="Catchability",main="Northern rockfish")
points(mcmc$M,mcmc$q1,col="yellow",pch="*",cex=1)

###### What else needs to be done?

