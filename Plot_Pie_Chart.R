#read the file that contains top 500 significant genes
raw_data <- read.table(file="clipboard",sep="\t",header=TRUE)
library(dplyr)
#Add a new column called p-adj
raw_data$p_adj <- p.adjust(raw_data[["pvalue"]], method = "fdr")
raw_data[,'Column1']=NULL
raw_data <- filter(raw_data, p_adj != "N/A")
#use FDR to select the threshold
df <- data.frame(matrix(NA,nrow=nrow(raw_data),ncol=ncol(raw_data)))

j=0

for(i in 1:nrow(raw_data)){
  #print(i)
  if(raw_data[i,12]<0.05){
    print(i)
    j=j+1
    df[j,]<-raw_data[i,]
  }else{
    
  }
}
df<-filter(df,X1!="N/A")
down_regulated_genes=1
up_regulated_genes=1

for(i in 1:nrow(df)){
  if(df[i,3]>0){
    up_regulated_genes=up_regulated_genes+1
  }else{
    down_regulated_genes=down_regulated_genes+1
  }
}

down_regulated_genes=down_regulated_genes-1
up_regulated_genes=up_regulated_genes-1

slices <- c(down_regulated_genes,up_regulated_genes)
lbls <- c("downregulated genes(269)","upregulated genes(546)")

pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="")
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Differentially Expressed Genes CCI vs SNI Using FDR 5%")

#for the significant inflammatory pathways
raw_data <- read.table(file="clipboard",sep="\t",header=TRUE)
raw_data$p_adj <- p.adjust(raw_data[["pval"]], method = "fdr")
raw_data[,'padj']=NULL
raw_data <- filter(raw_data, p_adj != "N/A")
a<-raw_data[grep("inflammatory", raw_data$desc), ]
b<-raw_data[grep("immune", raw_data$desc), ]
total<-rbind(a,b)

mice<-filter(raw_data,host=='Mmus')
mics<-mice[grep("inflammatory",mice$desc),]
mics<-filter(mics,MeSH=='A08.800.350.340')

rats<-filter(rats,desc=='inflammatory response')
time<-c(336,72,336,504,336,336,168)
rats[,'t']<-as.numeric(time)
