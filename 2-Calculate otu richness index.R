#### Jie Hu ####
#### this script will calculate diversity index for resident community 
## loading packages 
library(vegan) # only the ones that you need
library(BiodiversityR) 
library(lattice)
library(permute)
library(tcltk)
library(nlme)
library(picante)
library(ape)
### loading PAE function, PAE: phylogenetic-Abundance Evenness
PAE<-function(abund,tree){
  terms <- tree$edge[, 2] <= Ntip(tree)
  terminal.edges <- tree$edge.length[terms]
  names(terminal.edges)<-tree$tip.label
  PDvalue<-sum(tree$edge.length)
  spnames<-names(terminal.edges)
  termabun<-numeric(length(terminal.edges))
  for(i in 1:length(terminal.edges)){
    sp<-spnames[i] ## get species names
    termabun[i]<-(terminal.edges[sp])*(abund[sp]-1)
  }
  PAEdividend<-(sum(unlist(termabun)))+PDvalue
  PAEdivisor<-PDvalue+(((sum(abund[which(abund>=0)])/length(abund[which(abund>=0)]))-1)*sum(terminal.edges))
  PAEquotient<-PAEdividend/PAEdivisor
  return(PAEquotient)
}
## reading data 
 
tree1<-read.tree("otu_rep.2.tre",text=NULL,tree.names=NULL,skip=0,comment.char="#")

## subset data ##

ngs.otu <- as.data.frame(t(ngs.b.filter))      ## transposes
ngs.otu.01 <- as.data.frame(t(ngs.01))         ## transposes

#### calculate diversity index ####
otu.shannon<- diversity(ngs.otu, index = "shannon") 
otu.richness<-rowSums(ngs.otu>0)
otu.evenness<- otu.shannon/log(otu.richness)                  ## Pielou's evenness
otu.alpha<-fisher.alpha(ngs.b.filter,MARGIN=2)

n<-rowSums(ngs.otu==1)                                        ## number of singletons
m<-rowSums(ngs.otu==2)                                        ## number of doubletons
otu.Chao1<-otu.richness+n*(n-1)/(2*(m+1))                     ## calculate chao1 estimate 

#check that same row names are used in OTU table and matrix 
cophen.distb<-cophenetic(tree1)                               ## tree from company
tokeep.d<-rownames(ngs.b.filter) %in% rownames(cophen.distb)  ## basic knowledge
otu.nor.f<-ngs.b.filter[tokeep.d,]

tokeep.e<- rownames(cophen.distb) %in% rownames(ngs.b.filter) ## basic knowledge
cophen.distc<-cophen.distb[tokeep.e,tokeep.e]                 ## kind of distance matrix needed by the ses.mpd
otu.nor.f<-otu.nor.f[order(rownames(cophen.distc)),]

## calculate phylogenetic diversity
pd.1<-pd(t(otu.nor.f), tree1, include.root=F)
otu.pae<-t(otu.nor.f)

pae.list<-list()
for (i in 1:dim(otu.pae)[1]) pae.list[i]<-PAE(otu.pae[i,],tree1)

pae<-t(data.frame(pae.list))
  
#### calculate normalized otu richness  ####
otu.richness.ctr<-mean(otu.richness[env$div==0])
otu.richness.delta<-otu.richness-otu.richness.ctr
otu.richness.norm<-otu.richness.delta/otu.richness.ctr*100 
env$rel.ab<-10^env$phlD20/10^env$X16S ## calculate phlD gene relative abundance

newdata<-cbind(env,otu.richness,otu.shannon,pd.1,otu.alpha,otu.evenness,otu.richness.delta,otu.richness.norm)
newdata1<-newdata[1:48,]

#write.table(newdata,file="newdata.txt",sep="\t")

plot(newdata$otu.alpha~newdata$div)

