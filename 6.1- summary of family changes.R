#### To draa, at family level
ngs.family<-aggregate(ngs.b.filter, list(taxo.i$family), FUN = sum)     ## otu richness inside family, based on 0/1 matrix
ngs.family.relat<-aggregate(ngs.relat, list(taxo.i$family), FUN = sum)     ## otu richness inside family, based on 0/1 matrix
rownames(ngs.family.relat)<-ngs.family.relat[,1]
ngs.family.relat<-ngs.family.relat[order(rownames(ngs.family.relat)),]
ngs.family.relat.ck<-rowSums(ngs.family.relat[,2:49])/48

ngs.family.b <- ngs.family[, -1]
rownames(ngs.family.b)<-ngs.family[,1]        ## change column name
otu.in.family<-as.data.frame(t(ngs.family.b[,1:48]))
env<-env[1:48,]

#### calculate OTU abundance V.S. div and phlD at different level, based on normalized OTU table ####
#### analysis at family level
family.vs.div.Pvalue <- 0
family.vs.div.corvect <- 0

for (i in 1:dim(otu.in.family)[2])
{
  family.vs.div.Pvalue[i] <- (summary(lm(otu.in.family[,i]~env$div)))$coefficients[2, 4]
  family.vs.div.corvect[i] <- cor(otu.in.family[, i], env$div)
}
family.vs.div<-cbind(family.vs.div.Pvalue, family.vs.div.corvect)

rownames(family.vs.div) <- ngs.family[,1]
for (i in 1:dim(family.vs.div)[1]) 
  if (family.vs.div[i,1] == "NaN") family.vs.div[i,1] <-1
for (i in 1:dim(family.vs.div)[1]) 
  if (is.na(family.vs.div[i,2]) =="TRUE")  family.vs.div[i,2]<-0

family.color<-list()
for (i in 1:dim(family.vs.div)[1])
{
  if (family.vs.div[i,1]>= 0.05)
  {family.color[i]<- "grey"}
  else
  {family.color[i]<- "red"}
}

family.color<-as.data.frame(t(family.color))
family.color<-as.data.frame(t(family.color))
colnames(family.color)<-"color"
family.vs.div<-cbind(family.vs.div,family.color)

family.vs.div<-family.vs.div[order(family.vs.div[,2]),]

otu.in.family.rank<-as.data.frame(ngs.family.b)
otu.in.family.rank$Sum<-rowSums(ngs.family.b[,1:52])

otu.in.family.rank<-otu.in.family.rank[order(otu.in.family.rank[,5]),]

family.vs.div<-family.vs.div[order(rownames(family.vs.div)),]
otu.in.family.rank<-otu.in.family.rank[order(rownames(otu.in.family.rank)),]

family.rank<-cbind(otu.in.family.rank,family.vs.div,ngs.family.relat.ck)

family.rank<-family.rank[order(family.rank[,57]),]

rank <- seq(308,1)
family.rank<-cbind(family.rank,rank)

plot(family.rank$family.vs.div.corvect~rank)
library(ggplot2)

pdf("Figure 1F summary of family change.pdf", height=6, width=8)

family.change<-ggplot(data=family.rank,aes(x=log10(ngs.family.relat.ck), y=family.vs.div.corvect))+
  geom_point(size=3,color=family.rank$color)+
  geom_text(aes(label=ifelse(family.vs.div.Pvalue<0.05,row.names(family.rank),"")),hjust=0,vjust=0)+
  geom_hline(aes(yintercept=0),colour="#990000", linetype="dashed")+
  labs(x="Relative family abundance (log10 scale)", 
       y="Correlation coefficient of family abundance and the richness of introduced consortia")
#+theme_zg()

family.change
dev.off()

write.table(ngs.family.relat.ck, file="ngs.family.relat.ck.txt",row.names=T, dec=".", sep="\t")
