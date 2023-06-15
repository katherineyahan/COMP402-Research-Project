#Only keep the columns that we are interested in

#Eliminate the duplicate measurements
raw_data <- read.table(file="clipboard",sep="\t",header=TRUE)
library(dplyr)
raw_data <- distinct(raw_data,Gene, .keep_all=TRUE)
library("writexl")
write_xlsx(raw_data,"D:/COMP402 project2/DRG/Rat1_without_duplicates.xlsx")

#Go to gprofiler to humanize the genes
#Only keep the columns we want
humanized_data<- read.table(file="clipboard",sep="\t",header=TRUE)

#Eliminate the duplicates
humanized_data <- distinct(humanized_data,HumanGene, .keep_all=TRUE)
humanized_data <- distinct(humanized_data,Gene, .keep_all=TRUE)

#Eliminate the rows with N/A in HumanGene
humanized_data <- filter(humanized_data, HumanGene != "N/A")

#merge raw_data with humanized_data by Gene name

#If it is mouse, do this
raw_data[,"Gene"]=toupper(raw_data[,"Gene"])

dplyr::right_join(raw_data, humanized_data, by = "Gene", ignore_case=TRUE) -> Merged_table

#Drop the column that has the original gene names
Merged_table <-Merged_table[,-1]

#merge the file with Gene Length
GeneLength <- read.table(file="clipboard",sep="\t",header=TRUE)
dplyr::right_join(Merged_table, GeneLength, by = "HumanGene", ignore_case=TRUE) -> Merged_table
Merged_table <- filter(Merged_table, Sample1 != "N/A")

df <- data.frame(matrix(NA,nrow=nrow(Merged_table),ncol=(ncol(raw_data))))
col_names=list('Gene','Sample1','Sample2','Sample3','')
names(df)<-c(col_names)

#Change this accordingly
df[1]=Merged_table[,4]


#Convert from read counts to TPM
for (j in 1:(ncol(Merged_table)-2)){
  print(j)
  RPK<-list()
  for(i in 1:nrow(Merged_table)){
    
    each_RPK=Merged_table[i,j]/Merged_table[i,"length"]
    RPK<-append(RPK,each_RPK)
  }
  scaling_factor<-sum(unlist(RPK))/1000000
  TPM<-list()
  for(i in unlist(RPK)){
    each_TPM<-i/scaling_factor
    TPM<-append(TPM,each_TPM)
  }
  df[(j+1)] <- c(unlist(TPM))
}
write.table(df,"D:/COMP402 project2/SC/Mouse/TPM_GSE106803.txt",sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#Convert to FPKM
for(j in 1:ncol(raw_data)){
  print(j)
  #Calculate the per million scaling factor
  scaling_factor <- sum(Merged_table[,j])/1000000
  print(scaling_factor)
  FPKM <- list()
  #calculate the RPM for each count
  for(i in 1:nrow(Merged_table)){
    print(i)
    each_RPM=Merged_table[i,j]/scaling_factor
    each_FPKM=each_RPM/Merged_table[i,"length"]
    FPKM <-append(FPKM,each_FPKM)
  }
  df[(j+1)] <- c(unlist(FPKM))
}
write_xlsx(df,"C:/Users/yahan/Desktop/COMP402 project2/FPKM_Mouse2.xlsx")


fpkm_file<-read.table(file="clipboard",sep="\t",header=TRUE)
#Convert from RPKM/FPKM to TPM
df <- data.frame(matrix(NA,nrow=nrow(fpkm_file),ncol=(ncol(fpkm_file))))
col_names=list('Gene','Sample1','Sample2','Sample3')
names(df)<-c(col_names)

df[1]=fpkm_file[,1]
for(i in 1:(ncol(Merged_table)-1)){
  scaling_factor <- sum(Merged_table[,i])
  TPM<-list()
  for(j in 1:nrow(Merged_table)){
    each_TPM=(Merged_table[j,i]/scaling_factor)*100000
    TPM<-append(TPM,each_TPM)
  }
  df[(i+1)]<-c(unlist(TPM))
}
library("writexl")
write_xlsx(df,"D:/COMP402 project2/DRG/TPM_Mouse3.xlsx")
