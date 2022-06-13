library(tidyverse)
library(DESeq2)
#import data
setwd("D:/yyy/Merge_HUADA_YK_HUADA2/")
mycounts<-read.table("Mergedata_FOB.txt",header = TRUE,row.names = 1,sep ="\t" )
mycounts_new<-mycounts
a<-colnames(mycounts)
c<-data.frame(0,0,0,0,0,0)
colnames(c)<-c("j","k","p","q","y","z")

for (j in 1:9) {
  for (k in (j+1):10) {
    for (p in 11:26) {
      for (q in (p+1):27) {
        mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q])]
        condition<-factor(c(rep("KO",2),rep("WT",2)),levels = c("WT","KO"))
        colData<-data.frame(row.names = colnames(mycounts),condition)
        
        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
        dds <- DESeq(dds)
        
        res= as.data.frame(results(dds))
        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
        y<-length(rownames(res_hi))
        z<-length(rownames(res_lo))
        prdm1<-res[which(rownames(res) == "Prdm1"),]
        prdm1_a<-prdm1$log2FoldChange
        prdm1_b<-prdm1$padj
        if(y>z & (y/z)>=4 & prdm1_b<0.05 & prdm1_a>0.58 ){
          b<-as.data.frame(t(as.data.frame(c(j,k,p,q,y,z))))
          colnames(b)<-c("j","k","p","q","y","z")
          c<-rbind(c,b)
          write.table(mycounts,file = paste0("mycounts_YKHUADA",j,k,p,q,y,z,prdm1_a,".txt"),sep = "\t")
          write.table(res,file = paste0("mycounts_YKHUADA_Diff",j,k,p,q,y,z,prdm1_a,".txt"),sep = "\t")
        }      
        
      }
    }
    
    
  }
}
