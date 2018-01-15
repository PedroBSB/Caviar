RV1<- read.csv2("Data/RV.csv",header=F)
RV1<-RV1[1,c(-2,-3,-4,-5)]
RV2<- read.csv2("Data/RV.csv",header=F)
RV2<-RV2[2,c(-2,-3,-4,-5)]
RV<- read.csv2("Data/RV.csv",skip = 2,header=F)
RV<- RV[,c(-2,-3,-4,-5)]
low<-grepl('Trade Low',as.character(t(RV2)))
high<-grepl('Trade High',as.character(t(RV2)))

RV.low<-RV[,low]
RV.low<-apply(RV.low,2,function(x)as.numeric(as.character(x)))
RV.high<-RV[,high]
RV.high<-apply(RV.high,2,function(x)as.numeric(as.character(x)))
RV.final<- as.data.frame(RV.high-RV.low)
RV.final<-cbind(RV[,1],RV.final)
colnames(RV.final)<-c("Day",as.character(t(RV1[low])))
write.csv2(RV.final, "Data/indexRV.csv", row.names = F)
