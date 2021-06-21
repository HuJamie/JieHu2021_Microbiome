#### To draa, at genus level
ngs.genus<-aggregate(ngs.b.filter, list(taxo.i$genus), FUN = sum)     ## otu richness inside genus, based on 0/1 matrix
ngs.genus.relat<-aggregate(ngs.relat, list(taxo.i$genus), FUN = sum)     ## otu richness inside genus, based on 0/1 matrix
rownames(ngs.genus.relat)<-ngs.genus.relat[,1]
ngs.genus.relat<-ngs.genus.relat[order(rownames(ngs.genus.relat)),]
ngs.genus.relat.ck<-rowSums(ngs.genus.relat[,2:49])/48

ngs.genus.b <- ngs.genus[, -1]
rownames(ngs.genus.b)<-ngs.genus[,1]        ## change column name
otu.in.genus<-as.data.frame(t(ngs.genus.b[,1:48]))
env<-env[1:48,]

ngs.genus.b<-subset(ngs.genus.b,rownames(ngs.genus.b)!="Ralstonia")
ngs.rals$Ralstonia

#### calculate OTU abundance V.S. div and phlD at different level, based on normalized OTU table ####
#### analysis at genus level
genus.vs.div.Pvalue <- 0
genus.vs.div.corvect <- 0

for (i in 1:dim(otu.in.genus)[2])
{
  genus.vs.div.Pvalue[i] <- (summary(lm(otu.in.genus[,i]~env$div)))$coefficients[2, 4]
  genus.vs.div.corvect[i] <- cor(otu.in.genus[, i], env$div)
}
genus.vs.div<-cbind(genus.vs.div.Pvalue, genus.vs.div.corvect)

rownames(genus.vs.div) <- ngs.genus[,1]
for (i in 1:dim(genus.vs.div)[1]) 
  if (genus.vs.div[i,1] == "NaN") genus.vs.div[i,1] <-1
for (i in 1:dim(genus.vs.div)[1]) 
  if (is.na(genus.vs.div[i,2]) =="TRUE")  genus.vs.div[i,2]<-0

genus.color<-list()
for (i in 1:dim(genus.vs.div)[1])
{
  if (genus.vs.div[i,1]>= 0.05)
  {genus.color[i]<- "grey"}
  else
  {genus.color[i]<- "red"}
}

genus.color<-as.data.frame(t(genus.color))
genus.color<-as.data.frame(t(genus.color))
colnames(genus.color)<-"color"
genus.vs.div<-cbind(genus.vs.div,genus.color)

genus.vs.div<-genus.vs.div[order(genus.vs.div[,2]),]

otu.in.genus.rank<-as.data.frame(ngs.genus.b)
otu.in.genus.rank$Sum<-rowSums(ngs.genus.b[,1:52])

otu.in.genus.rank<-otu.in.genus.rank[order(otu.in.genus.rank[,5]),]

genus.vs.div<-genus.vs.div[order(rownames(genus.vs.div)),]
otu.in.genus.rank<-otu.in.genus.rank[order(rownames(otu.in.genus.rank)),]

genus.rank<-cbind(otu.in.genus.rank,genus.vs.div,ngs.genus.relat.ck)

genus.rank<-genus.rank[order(genus.rank[,57]),]

rank <- seq(581,1)
genus.rank<-cbind(genus.rank,rank)

plot(genus.rank$genus.vs.div.corvect~rank)
library(ggplot2)

pdf("Figure 1F summary of genus change.pdf", height=10, width=12)

genus.change<-ggplot(data=genus.rank,aes(x=log10(ngs.genus.relat.ck), y=genus.vs.div.corvect))+
  geom_point(size=3,color=genus.rank$color)+
  geom_text(aes(label=ifelse(genus.vs.div.Pvalue<0.05,row.names(genus.rank),"")),hjust=0,vjust=0)+
  geom_hline(aes(yintercept=0),colour="#990000", linetype="dashed")+
  labs(x="Relative genus abundance (log10 scale)", 
       y="Correlation coefficient of genus abundance and the richness of introduced consortia")
#+theme_zg()

genus.change
dev.off()

write.table(ngs.genus.relat.ck, file="ngs.genus.relat.ck.txt",row.names=T, dec=".", sep="\t")
