#Read the file
human1<-read.table(file="clipboard",sep="\t",header=TRUE)
human2<-read.table(file="clipboard",sep="\t",header=TRUE)
human3<-read.table(file="clipboard",sep="\t",header=TRUE)
human4<-read.table(file="clipboard",sep="\t",header=TRUE)

mouse1<-read.table(file="clipboard",sep="\t",header=TRUE)
mouse2<-read.table(file="clipboard",sep="\t",header=TRUE)
mouse4<-read.table(file="clipboard",sep="\t",header=TRUE)

rat1<-read.table(file="clipboard",sep="\t",header=TRUE)
rat2<-read.table(file="clipboard",sep="\t",header=TRUE)

#Merge them
library(dplyr)
dplyr::right_join(human1, human2, by = "Gene", ignore_case=TRUE) -> Merged_Human
dplyr::right_join(Merged_Human,human3, by="Gene", ignore_case=TRUE)->Merged_Human
dplyr::right_join(Merged_Human,human4, by="Gene", ignore_case=TRUE)->Merged_Human
col_names=list('Gene','Sample1','Sample2','Sample3','Sample4','Sample5','Sample6','Sample7','Sample8','Sample9','Sample10','Sample11','Sample12','Sample13')
names(Merged_Human)<-c(col_names)
Merged_Human <- filter(Merged_Human, Sample1.x != "N/A")

dplyr::right_join(mouse1, mouse2, by = "Gene", ignore_case=TRUE) -> Merged_Mouse
dplyr::right_join(Merged_Mouse, mouse4, by = "Gene", ignore_case=TRUE) -> Merged_Mouse
col_names=list('Gene','Sample1','Sample2','Sample3','Sample4','Sample5','Sample6','Sample7','Sample8','Sample9')
names(Merged_Mouse)<-c(col_names)
Merged_Mouse <- filter(Merged_Mouse, Sample1 != "N/A")

dplyr::right_join(rat1, rat2, by = "Gene", ignore_case=TRUE) -> Merged_Rat
col_names=list('Gene','Sample1','Sample2','Sample3','Sample4','Sample5','Sample6','Sample7','Sample8','Sample9')
names(Merged_Rat)<-c(col_names)
Merged_Rat <- filter(Merged_Rat, Sample1 != "N/A")

#Calculate the correlation and plot the distribution
library(MASS)
Merged_Human[,45] <- NULL
HH<-list()
sample_Human<-list()
for (i in 2:(ncol(Merged_Human)-1)){
  for(j in (i+1):ncol(Merged_Human)){
    #print(i)
    #print(j)
    cor<-cor(Merged_Human[,i],Merged_Human[,j],method="spearman",use="complete.obs")
    HH<-append(HH,cor)
    sample<-c(i,j)
    sample_Human<-append(sample_Human,sample)
    if (cor < 0.5){
      print(i)
      print(j)
    }
  }
}
densHH<-density(unlist(HH))
plot(densHH,xlab='spearman correlation',main='Distribution of Human and Human Correlation')
abline(v=mean(unlist(HH)),col="red")
#text("Mean",x=v,y=max(y),srt=-90,pos=4)
legend("topleft",legend="Mean",pch="|",col="red")
abline(v=median(unlist(HH)),col="blue")
#text("Median",x=v,y=max(y),srt=-90,pos=4)
legend("topright",legend="Median",pch="|",col="blue")

MM<-list()
sample_Mouse<-list()
for (i in 2:(ncol(Merged_Mouse)-1)){
  for(j in (i+1):ncol(Merged_Mouse)){
    cor<-cor(Merged_Mouse[,i],Merged_Mouse[,j],method="spearman",use="complete.obs")
    MM<-append(MM,cor)
    sample<-c(i,j)
    sample_Mouse<-append(sample_Mouse,sample)
  }
}
densMM<-density(unlist(MM))
plot(densMM,xlab='spearman correlations',main='Distribution of Mouse and Mouse Correlation')
abline(v=mean(unlist(MM)),col="red")
#text("Mean",x=v,y=max(y),srt=-90,pos=4)
legend("topleft",legend="Mean",pch="|",col="red")
abline(v=median(unlist(MM)),col="blue")
#text("Median",x=v,y=max(y),srt=-90,pos=4)
legend("topright",legend="Median",pch="|",col="blue")


RR<-list()
sample_Rat<-list()
for (i in 2:(ncol(Merged_Rat)-1)){
  for(j in (i+1):ncol(Merged_Rat)){
    cor<-cor(Merged_Rat[,i],Merged_Rat[,j],method="spearman",use="complete.obs")
    RR<-append(RR,cor)
    sample<-c(i,j)
    sample_Rat<-append(sample_Rat,sample)
  }
}
densRR<-density(unlist(RR))
plot(densRR,xlab='spearman correlations',main='Distribution of Rat and Rat Correlations')
abline(v=mean(unlist(RR)),col="red")
#text("Mean",x=v,y=max(y),srt=-90,pos=4)
legend("topleft",legend="Mean",pch="|",col="red")
abline(v=median(unlist(RR)),col="blue")
#text("Median",x=v,y=max(y),srt=-90,pos=4)
legend("topright",legend="Median",pch="|",col="blue")

merged <- dplyr::right_join(Merged_Human, Merged_Mouse, by = "Gene", ignore_case=TRUE)
#col=c("Gene","H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","M1","M2","R1","R2","R3","R4","R5")
#colnames(merged)=col
merged <- merged[!(merged$Sample1.x=="N/A"),]

library(dplyr)

Hsap <- merged[,2:69]
Mmus <- merged[,70:78]
HM<-list()
for (i in 1:ncol(Hsap)){
  print(i)
  for(j in 1:ncol(Mmus)){
    print(j)
    corr=cor(Hsap[,i],Mmus[,j],method="spearman",use="complete.obs")
    HM<-append(HM,corr)
  }
}
densHM<-density(unlist(HM))
plot(densHM,xlab='spearman correlations',main='Distribution of Human and Mouse Correlations')
abline(v=mean(unlist(HM)),col="red")
legend("topleft",legend="Mean",pch="|",col="red")
abline(v=median(unlist(HM)),col="blue")
legend("topright",legend="Median",pch="|",col="blue")

merged <- dplyr::right_join(Merged_Human, Merged_Rat, by = "Gene", ignore_case=TRUE)
merged <- merged[!(merged$Sample1.x=="N/A"),]

library(dplyr)

Hsap <- merged[,2:69]
Rnor <- merged[,70:78]
HR<-list()
for (i in 1:ncol(Hsap)){
  print(i)
  for(j in 1:ncol(Rnor)){
    print(j)
    corr=cor(Hsap[,i],Rnor[,j],method="spearman",use="complete.obs")
    HR<-append(HR,corr)
  }
}
densHR<-density(unlist(HR))
plot(densHR,xlab='spearman correlations',main='Distribution of Human and Rat Correlations')
abline(v=mean(unlist(HR)),col="red")
legend("topleft",legend="Mean",pch="|",col="red")
abline(v=median(unlist(HR)),col="blue")
legend("topright",legend="Median",pch="|",col="blue")

merged <- dplyr::right_join(Merged_Mouse, Merged_Rat, by = "Gene", ignore_case=TRUE)
merged <- merged[!(merged$Sample1.x=="N/A"),]

library(dplyr)

Rnor <- merged[,9:17]
Mmus<-merged[,2:8]
MR<-list()
for (i in 1:ncol(Mmus)){
  print(i)
  for(j in 1:ncol(Rnor)){
    print(j)
    corr=cor(Mmus[,i],Rnor[,j],method="spearman",use="complete.obs")
    MR<-append(MR,corr)
  }
}
densMR<-density(unlist(MR))
plot(densMR,xlab='spearman correlations',main='Distribution of Mouse and Rat Correlations')
abline(v=mean(unlist(MR)),col="red")
legend("topleft",legend="Mean",pch="|",col="red")
abline(v=median(unlist(MR)),col="blue")
legend("topright",legend="Median",pch="|",col="blue")

library('plotrix')
z=(abs(mean(unlist(HM))-mean(unlist(HR))))/sqrt(std.error(unlist(HM))^2+std.error(unlist(HR))^2)
p=2*pnorm(-abs(z))
