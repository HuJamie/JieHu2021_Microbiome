####
library(ggplot2)
library(grid)
library(gridExtra)
#### set up backgroud as white for ggplot2 ####
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black', vjust=0.1),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))
}

#######
pdf("Figure S1 Total bacterial density .pdf", height=3, width=4)
fig.16s<-ggplot(newdata,aes(x=as.factor(div), y=X16S))+geom_boxplot()+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader richness", y="Total bacterial density in the resident community")+theme_zg()
fig.16s
dev.off()

lm.s2<-aov(X16S~as.factor(div),newdata)
summary(lm.s2)
TukeyHSD(lm.s2)
plot(TukeyHSD(lm.s2))


pdf("Figure S2 Pseudomonas relative abundance .pdf", height=3, width=4)
fig.por<-ggplot(newdata,aes(x=as.factor(div), y=rel.ab*100))+geom_boxplot()+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader richness", y="Invader proportion (phlD gene copies/16S gene copies %)")+
  theme_zg()
fig.por
dev.off()

lm.s1<-aov(rel.ab*100~as.factor(div),newdata)
summary(lm.s1)
TukeyHSD(lm.s1)
plot(TukeyHSD(lm.s1))

############ How many Pseudomonas sequences inside amplicon library,by Hu Jie 
#### Reading data 
library(ggplot2)
library(grid)
####

sequ<-read.table("sequence99.txt", header=T, row.names=1, sep="\t")  ## sequence rawdata
env<-read.table("env.txt", sep="\t", header=T, row.names=1)        ## env

env.t<-env[1:48,5:12]
sequ1<-sequ[1:48,]

sequ.48<-env.t*sequ1
sequ.48$tol<-rowSums(sequ.48)
boxplot(sequ.48$tol~env[1:48,]$div)

sequ$tol<-rowSums(sequ)
boxplot(sequ$tol~env$div)
anova(lm(sequ.48$tol~env[1:48,]$div))

env1<-cbind(sequ,env)

pdf("Figure S1.B Pseudomonas relative abundance based on sequence.pdf", height=3, width=4)
fig.seq<-ggplot(env1,aes(x=as.factor(div), y=tol))+geom_boxplot()+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader richness", y="Relative sequence abundance")+theme_zg()
fig.seq
dev.off()

env.2<-cbind(env,sequ)
lm.div.seq<-lm(tol~div,env.2[1:48,])
anova(lm.div.seq)
summary(lm.div.seq)
lm.div.seq1<-aov(tol~as.factor(div),env.2[1:48,])
TukeyHSD(lm.div.seq1)


boxplot(sequ$MVP1.4~env$div)
boxplot(sequ$Q2.87~env$div)
boxplot(sequ$CHA0~env$div)
boxplot(sequ$F113~env$div)
boxplot(sequ$Phl1c2~env$div)
boxplot(sequ$Pf.5~env$div)
boxplot(sequ$M1.96~env$div)
boxplot(sequ$Q8R1.96~env$div)

sequ.t<-data.frame(t(sequ))

boxplot(sequ.t[1:8,]$s01~c(1,2,3,4,5,6,7,8))
boxplot(sequ.t[1:8,]$s45~c(1,2,3,4,5,6,7,8))
boxplot(sequ.t[1:8,]$s46~c(1,2,3,4,5,6,7,8))
boxplot(sequ.t[1:8,]$s02~c(1,2,3,4,5,6,7,8))

