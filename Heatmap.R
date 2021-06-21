#####heatmap od changed OTUs
#### To draa, at family level
ngs.relat<-ngs.relat[order(rownames(ngs.relat)),]
ngs.otu<-ngs.b.filter[order(rownames(ngs.b.filter)),]
ngs.relat.ck<-rowSums(ngs.relat[,1:48])/48

otu.ngs<-as.data.frame(t(ngs.otu[,1:48]))
env<-env[1:48,]

#### calculate OTU abundance V.S. div and phlD at different level, based on normalized OTU table ####
#### analysis at family level
OTU.vs.div.Pvalue <- 0
OTU.vs.div.corvect <- 0

for (i in 1:dim(otu.ngs)[2])
{
  OTU.vs.div.Pvalue[i] <- (summary(lm(otu.ngs[,i]~env$div)))$coefficients[2, 4]
  OTU.vs.div.corvect[i] <- cor(otu.ngs[, i], env$div)
}
otu.vs.div<-cbind(OTU.vs.div.Pvalue, OTU.vs.div.corvect)

rownames(otu.vs.div) <- rownames(ngs.otu)
for (i in 1:dim(otu.vs.div)[1]) 
  if (otu.vs.div[i,1] == "NaN") otu.vs.div[i,1] <-1
for (i in 1:dim(otu.vs.div)[1]) 
  if (is.na(otu.vs.div[i,2]) =="TRUE")  otu.vs.div[i,2]<-0

otu.color<-list()
for (i in 1:dim(otu.vs.div)[1])
{
  if (otu.vs.div[i,1]>= 0.05)
  {otu.color[i]<- "grey"}
  else
  {otu.color[i]<- "red"}
}

otu.color<-as.data.frame(t(otu.color))
otu.color<-as.data.frame(t(otu.color))
colnames(otu.color)<-"color"
otu.vs.div<-cbind(otu.vs.div,otu.color)

otu.vs.div<-otu.vs.div[order(otu.vs.div[,2]),]

otu.rank<-as.data.frame(ngs.otu)
otu.rank$Sum<-rowSums(ngs.otu[,49:52])

#otu.rank<-otu.rank[order(otu.rank[,53]),]

otu.vs.div<-otu.vs.div[order(rownames(otu.vs.div)),]
otu.rank<-otu.rank[order(rownames(otu.rank)),]

otu.rank<-cbind(otu.rank,otu.vs.div,ngs.relat.ck)

otu.rank<-otu.rank[order(otu.rank[,57]),]

rank <- seq(8631,1)
otu.rank<-cbind(otu.rank,rank)

otu.heat<-otu.rank[,54:58]
otu.heat<-subset(otu.heat,otu.heat$color=="red")
otu.heat.d<-subset(otu.heat,otu.heat$OTU.vs.div.corvect<0)
otu.heat.i<-subset(otu.heat,otu.heat$OTU.vs.div.corvect>0)

otu.d<-as.data.frame(otu.heat.d$ngs.relat.ck)
rownames(otu.d)<-rownames(otu.heat.d)
tokeep.d<-rownames(ngs.otu) %in% rownames(otu.d)  ## basic knowledge
heat.d<-ngs.otu[tokeep.d,1:48]

otu.i<-as.data.frame(otu.heat.i$ngs.relat.ck)
rownames(otu.i)<-rownames(otu.heat.i)
tokeep.i<-rownames(ngs.otu) %in% rownames(otu.i)  ## basic knowledge
heat.i<-ngs.otu[tokeep.i,1:48]
colnames(heat.i)<-env$div
heat.i.t<-as.data.frame(t(heat.i))
heat.i.t<-aggregate(heat.i.t, list(env$div), FUN = sum)     ## otu richness inside family, based on 0/1 matrix

colnames(heat.d)<-env$div
heat.d.t<-as.data.frame(t(heat.d))
heat.d.t<-aggregate(heat.d.t, list(env$div), FUN = sum)     ## otu richness inside family, based on 0/1 matrix


heat.i.t<-as.matrix(heat.i.t[,-1])
heat.d.t<-as.matrix(heat.d.t[,-1])

heatmap(log10(heat.i.t+1), Colv = NA, Rowv = NA, scale="column")
heatmap(log10(heat.d.t+1), Colv = NA, Rowv = NA, scale="column")

install.packages("gplots")
library("gplots")
heatmap.2(log2(heat.i.t+1),Colv=NA,Rowv=NA,scale="none",col=bluered(100), 
          trace="none",density.info="none")


