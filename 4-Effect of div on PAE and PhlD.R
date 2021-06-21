#####
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

####reading the data ####
##### Figure 3 effect of div on PD and PAE ####
mean.pd<-mean(newdata[49:52,]$PD)
mean.pae<-mean(newdata[49:52,]$pae)

pdf("Figure 1B effect of div on pae.pdf", height=3, width=4)

fig.pae.div<-ggplot(data=newdata1,aes(x=log2(div), y=pae))+geom_smooth(method="lm",span=1,se=FALSE)+
  geom_jitter(position = position_jitter(width =.2))+
  geom_hline(aes(yintercept=mean.pae),colour="#990000", linetype="dashed")+
  labs(x="Invader richness", y="Resident community phylogenetic abundance evenness")+
  theme_zg()
fig.pae.div

dev.off()

lm<-lm(PD~poly(div,2),newdata1)
anova(lm)
summary(lm)

lm.pae.div<-lm(pae~div,newdata1)
anova(lm.pae.div)
summary(lm.pae.div)

lm.pae.div1<-aov(pae~as.factor(div),newdata1)
TukeyHSD(lm.pae.div1)
