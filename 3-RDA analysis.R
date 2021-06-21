#### Jie Hu, Alex
#### This script will make RDA analyses from the community structure ###
#### libraries ####
library(BiodiversityR) # only the ones that you need
library(vegan)   

#### RDA of div+phlD ####
ngs.otu.01<-t(ngs.01)
rda1<-rda(ngs.otu.01[1:48,]~div+phlD20,env[1:48,]) 
rda1
summary(rda1,display="sites")

## important result ##

adonis.div1<-cca(t(ngs.01)~breadth+Gas+toxin1+toxin2+siderophore.SU.+IAA+p.solu+DAPG.ppm,env,by=NULL, method="euclidean", binary=T) 
adonis.div1

adonis.2<-adonis(t(ngs.b.filter)~toxin1+toxin2+breadth+Gas+siderophore.SU.+IAA+p.solu+DAPG.ppm,env,by=NULL, method="euclidean", binary=T) 
adonis.2

adonis.div2<-adonis(t(ngs.01)~phlD20+as.factor(div),env,by=NULL, method="euclidean", binary=T) 
adonis.div2

adonis.identity<-adonis(t(ngs.01[,1:48])~MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96,newdata1,by=NULL,method="euclidean", binary=T) 
adonis.identity

pdf("RDA of div and phlD.pdf",height=5,width=6)

palette1<-c("red1","red2","red3","red4","black")
plot.otu.rda <- plot(rda1, choices = c(1, 2), col="white",scaling = 2)
ordisymbol(plot.otu.rda,env,"div",col=1,colors=palette1, pchs=F, legend=T, legend.x="topright", legend.ncol=1)
ordiellipse(plot.otu.rda,as.factor(env$div),kind ="ehull",col=palette1, legend=F)

dev.off()

