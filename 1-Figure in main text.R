#############################
##### Jie Hu, 20180211, figures in main text
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

#### reading the data ####
#### What's happening ####
mean.pae<-mean(newdata[49:52,]$pae)
mean.nti<-mean(newdata[49:52,]$mntd.obs.z)

pdf("Figure 1F effect of div on PAE and multifunctionality.pdf", height=5, width=13)

fig.pae.div<-ggplot(data=newdata1,aes(x=log2(div), y=pae))+
  #geom_boxplot(aes(factor=factor(level)))+
  geom_smooth(method="lm",span=1,se=FALSE)+
  geom_jitter()+
  geom_hline(aes(yintercept=mean.pae),colour="#990000", linetype="dashed")+
  labs(x="Invader richness", y="Resident community phylogenetic abundance evenness")+
  theme_zg()

fig.div.multi<-ggplot(data=newdata.1,aes(x=log2(div), y=multi2))+
  stat_smooth(method="lm",fullrange = TRUE,se=FALSE)+geom_point(size=3,aes(color=factor(div)))+
  labs(x="Invader richnesss", y="Multifunctionality Index")+
  theme_zg()

grid.arrange(fig.pae.div,fig.div.multi,ncol=2)
dev.off()

pdf("Figure S1 effect of PAE on multifunctionality.pdf", height=5, width=7)

fig.pae.multi<-ggplot(data=newdata.1,aes(x=pae, y=multi2))+
  scale_x_continuous(limits=c(0.46,0.66)) +
  stat_smooth(method="lm",fullrange = TRUE,se=FALSE)+geom_point(size=3,aes(color=factor(div)))+
  labs(x="Resident community phylogenetic abundance evenness", y="Multifunctionality Index")+
  theme_zg()
fig.pae.multi
dev.off()


t.test(pae~as.factor(level),newdata1)
lm.pae.div<-aov(pae~as.factor(level),newdata1)
anova(lm.pae.div)
summary(lm.pae.div)
TukeyHSD(lm.pae.div)

lm.div.pae<-lm(pae~div,newdata.1)
anova(lm.div.pae)
summary(lm.div.pae)

lm.div.multi<-lm(multi2~log2(div),newdata.1)
anova(lm.div.multi)
summary(lm.div.multi)

lm.pae.multi<-lm(multi2~pae,newdata.1)
anova(lm.pae.multi)
summary(lm.pae.multi)

