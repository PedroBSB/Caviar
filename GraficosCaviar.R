setwd("C:\\Users\\Pedro Albuquerque\\Desktop\\R")
dados<-read.csv("hits_consolidado.csv")
dados<-dados[,-1]
dados$Model<-as.character(dados$Model)

dados$Model[dados$Model=='Adaptive']<-' Adaptive'
dados$Model[dados$Model=='symmetricAbs']<-' Symmetric Absolute Value'
dados$Model[dados$Model=='linearGARCH']<-' Symmetric Absolute Value'
dados$Model[dados$Model=='linearGARCHmu']<-' Symmetric Absolute Value with mu'
dados$Model[dados$Model=='linearTGARCH']<-' Asymmetric Slope'
dados$Model[dados$Model=='linearTGARCHmu']<-' Asymmetric Slope with mu'
dados$Model[dados$Model=='IndirectGARCH']<-' Indirect GARCH'
dados$Model[dados$Model=='standardGARCHmu']<-' Indirect GARCH with mu'
dados$Model[dados$Model=='GJRGARCH']<-' GJR-GARCH '
dados$Model[dados$Model=='GJRGARCHmu']<-' GJR-GARCH with mu'
dados$Model[dados$Model=='TCAV']<-' Threshold CAViaR'
dados$Model[dados$Model=='rangeValue']<-' Range Value '
dados$Model[dados$Model=='trigGARCH']<-' Threshold Range Indirect GARCH'
dados$Model[dados$Model=='trvGARCH']<-' Threshold Range Value'


train<-subset(dados,base=="train")
valid<-subset(dados,base=="valid")

train$lb<-train$Mean.Hits-2*sqrt(train$Var.Hits)
train$ub<-train$Mean.Hits+2*sqrt(train$Var.Hits)


valid$lb<-valid$Mean.Hits-2*sqrt(valid$Var.Hits)
valid$ub<-valid$Mean.Hits+2*sqrt(valid$Var.Hits)

library(ggplot2)

ggplot(data = train, aes(Model, Mean.Hits, color = Model)) +
  geom_point(size = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(
    aes(ymin = lb, ymax = ub),
    width = 0.1,
    linetype = "solid",
    position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept = 0.05), linetype = "dashed", color="red") + 
  theme_bw()+facet_wrap(~ativo,nrow=3)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("train.pdf", height=20, width=29, units='cm')



ggplot(data = valid, aes(Model, Mean.Hits, color = Model)) +
  geom_point(size = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(
    aes(ymin = lb, ymax = ub),
    width = 0.1,
    linetype = "solid",
    position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept = 0.05), linetype = "dashed", color="red") + 
  theme_bw()+facet_wrap(~ativo,nrow=3)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("valid.pdf", height=20, width=29, units='cm')

