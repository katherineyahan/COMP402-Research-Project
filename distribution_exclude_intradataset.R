#change the column names
colnames(df) <- c(new_col1_name,new_col2_name,new_col3_name)
#for loop to add the dataset number for the column names
for(i in 2:ncol(rat2)){
  colname=colnames(rat2)[i]
  new_name=paste("rat2_",colname)
  colnames(rat2)[i]<-new_name
}

#merged the files

#calculate the correlations and plot the distribution exclude intra-dataset correlations
MM<-list()
sample_Mouse<-list()
for (i in 2:(ncol(Merged_Mouse)-1)){
  for(j in (i+1):ncol(Merged_Mouse)){
    #extract the column names and compare->if does not have the same word before "-",calculate
    col1=colnames(Merged_Mouse)[i]
    col2=colnames(Merged_Mouse)[j]
    result1 <- strsplit(col1, "_")[[1]][1]
    result2 <- strsplit(col2, "_")[[1]][1]
    if(result1 != result2){
    cor<-cor(Merged_Mouse[,i],Merged_Mouse[,j],method="spearman",use="complete.obs")
    MM<-append(MM,cor)
    sample<-c(i,j)
    sample_Mouse<-append(sample_Mouse,sample)
    }
  }
}


RR<-list()
sample_Rat<-list()
for (i in 2:(ncol(Merged_Rat)-1)){
  for(j in (i+1):ncol(Merged_Rat)){
    col1=colnames(Merged_Rat)[i]
    col2=colnames(Merged_Rat)[j]
    result1 <- strsplit(col1, "_")[[1]][1]
    result2 <- strsplit(col2, "_")[[1]][1]
    if(result1 != result2){
    cor<-cor(Merged_Rat[,i],Merged_Rat[,j],method="spearman",use="complete.obs")
    RR<-append(RR,cor)
    sample<-c(i,j)
    sample_Rat<-append(sample_Rat,sample)
    }
  }
}


HH<-list()
sample_Human<-list()
for (i in 2:(ncol(Merged_Human)-1)){
  for(j in (i+1):ncol(Merged_Human)){
    col1=colnames(Merged_Human)[i]
    col2=colnames(Merged_Human)[j]
    result1 <- strsplit(col1, "_")[[1]][1]
    result2 <- strsplit(col2, "_")[[1]][1]
    if(result1 != result2){
      cor<-cor(Merged_Human[,i],Merged_Human[,j],method="spearman",use="complete.obs")
      HH<-append(HH,cor)
      sample<-c(i,j)
      sample_Human<-append(sample_Human,sample)
    }
  }
}