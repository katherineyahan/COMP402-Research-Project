#Calculate the pairwise correlation

#Read the TPM file
human <- read.table(file="clipboard",sep="\t",header=TRUE)
rat <- read.table(file="clipboard",sep="\t",header=TRUE)
mouse <- read.table(file="clipboard",sep="\t",header=TRUE)

#Only keep the gene and mean columns in each file
col<-list("Gene","Mean")
human<-data.frame(human$Gene,human$Mean)
names(human)<-c(col)
rat<-data.frame(rat$Gene,rat$Mean)
names(rat)<-c(col)
mouse<-data.frame(mouse$Gene,mouse$Mean)
names(mouse)<-c(col)

#merged them pairwise
dplyr::right_join(human, rat, by = "Gene", ignore_case=TRUE) -> H_R
dplyr::right_join(human, mouse, by = "Gene", ignore_case=TRUE) -> H_M
dplyr::right_join(rat, mouse, by = "Gene", ignore_case=TRUE) -> R_M

col<-list("Gene","Human","Rat")
names(H_R)<-c(col)

col<-list("Gene","Human","Mouse")
names(H_M)<-c(col)

col<-list("Gene","Rat","Mouse")
names(R_M)<-c(col)

#Calculate the correlation
H_R <- H_R[!(H_R$Human=="N/A"),]
H_R$Human->Hsap
H_R$Rat->Rnor
cor(Hsap,Rnor,method="spearman",use="complete.obs")

H_M <- H_M[!(H_M$Human=="N/A"),]
H_M$Human->Hsap
H_M$Mouse->Mmus
cor(Hsap,Mmus,method="spearman",use="complete.obs")

R_M <- R_M[!(R_M$Rat=="N/A"),]
R_M$Mouse->Mmus
R_M$Rat->Rnor
cor(Mmus,Rnor,method="spearman",use="complete.obs")

#Select only the protein encoding genes
genes_group <- read.table(file="clipboard",sep="\t",header=TRUE)
columns = c("Genes","Human","Rat")
df <- data.frame(matrix(ncol=length(columns),nrow=nrow(H_R)))
colnames(df)=columns
df2<-data.frame(matrix(ncol=length(columns),nrow=nrow(H_R)))
colnames(df2)=columns
protein_genes <- list()
for(j in 1:nrow(genes_group)){
  if (genes_group[j,2]=="protein coding"){
    protein_genes<-append(protein_genes,genes_group[j,1])
  }
}

row1 <- 1
row2<-1
for (i in 1:nrow(H_R)){
  if (H_R[i,1] %in% protein_genes){
    df[row1,1]<-H_R[i,1]
    df[row1,2]<-H_R[i,2]
    df[row1,3]<-H_R[i,3]
    row1<-row1+1
  }else{
    df2[row2,1]<-H_R[i,1]
    df2[row2,2]<-H_R[i,2]
    df2[row2,3]<-H_R[i,3]
    row2<-row2+1
  }
}

#Correlation of protein encoding genes
df$Human->Hsap
df$Rat->Rnor
df$Mouse->Mmus
cor(Hsap,Mmus,method="spearman",use="complete.obs")
cor(Hsap,Rnor,method="spearman",use="complete.obs")
cor(Rnor,Mmus,method="spearman",use="complete.obs")

cor(Hsap,Mmus,method="pearson",use="complete.obs")
cor(Hsap,Rnor,method="pearson",use="complete.obs")

#Correlation of non protein genes
df2$Human->Hsap
df2$Mouse->Mmus
df2$Rat->Rnor
cor(Hsap,Mmus,method="spearman",use="complete.obs")
cor(Hsap,Rnor,method="spearman",use="complete.obs")
cor(Rnor,Mmus,method="spearman", use="complete.obs")

cor(Hsap,Mmus,method="pearson",use="complete.obs")
cor(Hsap,Rnor,method="pearson",use="complete.obs")

