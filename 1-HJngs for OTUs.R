############
############ Diverse invaders manuscript script,by Hu Jie 
#### Reading data 
ngs.r<-read.table("otu.raw.txt", header=T, row.names=1, sep="\t")  ## otu rawdata
env<-read.table("env.txt", sep="\t", header=T, row.names=1, dec = ".")        ## environmental properties

ngs.i<-ngs.r[,1:dim(env)[1]]                                   ## only abundance
ngs.i<-ngs.i[,order(colnames(ngs.i))]                          ## order the data 
taxo.i<-ngs.r[,-c(1:dim(env)[1])]                              ## only taxonomy

ngs.8<-ngs.r[,45:48]  
ngs.8$sum<-rowSums(ngs.8)

ngs.8<-ngs.8[order(ngs.8$sum),]
write.table(ngs.8, file="ngs.8.txt",row.names=T, dec=".", sep="\t")

#### calculate value for method description 
colSums(ngs.i)
mean(colSums(ngs.i))
min(colSums(ngs.i))
max(colSums(ngs.i))

ngs.ra<-ngs.i[10:15,]
rarecurve(ngs.ra,step=1,sample=28809, xlab="Sample Size", ylab="Species",col="blue",label=FALSE)

env$nor.16<-10^(env$X16S)*0.3/50
env$nor.16<-10^(env$X16S)

## calculate CV for sequencing depth and qPCR of 16S data
mean.16s<-mean(env$nor.16)
sd.16s<-sd(env$nor.16)
cv.16s<-sd.16s/mean.16s

mean.seq<-mean(colSums(ngs.i))
sd.seq<-sd(colSums(ngs.i))
cv.seq<-sd.seq/mean.seq

mean.16s.f<-mean(colSums(ngs.b.filter))
min(colSums(ngs.b.filter))
max(colSums(ngs.b.filter))
sd.16s.f<-sd(colSums(ngs.b.filter))
cv.16s.f<-sd.16s.f/mean.16s.f

mean.seq.f<-mean(colSums(ngs.b.filter1))
min(colSums(ngs.b.filter1))
max(colSums(ngs.b.filter1))
sd.seq.f<-sd(colSums(ngs.b.filter1))
cv.seq.f<-sd.seq.f/mean.seq.f

#mean.gr<-aggregate(env$nor.16, FUN=mean,by=list(env$plant))
#sd.gr<-aggregate(grow.24$TSA0.1.od, FUN=sd,by=list(grow.24$plant))

#### calculate relative abundance of each otu 
ngs.relat<-0*ngs.i                                                       ## creat a empty matrix
for(i in 1:dim(ngs.relat)[2]) ngs.relat[,i]<-ngs.i[,i]/env$nor.16[i] ## loop to calculate relative abundance
colSums(ngs.relat)                                                       ## check if calculate correctly

#### calculate relative abundance of each otu 
ngs.relat<-0*ngs.i                                                       ## creat a empty matrix
for(i in 1:dim(ngs.relat)[2]) ngs.relat[,i]<-ngs.i[,i]/colSums(ngs.i)[i] ## loop to calculate relative abundance
colSums(ngs.relat)                                                       ## check if calculate correctly

#### Standardization of abundance
ngs.b<-round(ngs.relat*min(colSums(ngs.i)),0) ## recover otu number based relative otu abundance
#ngs.b<-round(ngs.relat*(env$nor.16),0) ## recover otu number based relative otu abundance
ngs.rt<-data.frame(ngs.b,taxo.i)              ## combine otu and taxonomy
ngs.b.filter<-ngs.rt[,1:dim(env)[1]]          ## subset based on dimensionality of env.txt

ngs.b$zero.num<-rowSums(ngs.b==0)             ## calculate the zero amount of each row
ngs.b$sum<-rowSums(ngs.b)                     ## calculate sum of each row, it makes some mistakes!!!!

colb<-data.frame(colSums(ngs.b!=0))           ## OTU richness of each sample

ngs.01<-as.data.frame(1*(ngs.b.filter>0))     ## transfer to 0/1 matrix 

