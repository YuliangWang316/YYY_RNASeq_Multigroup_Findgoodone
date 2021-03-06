library(tidyverse)
library(DESeq2)
#import data
setwd("E:/YYY_HUADA/counts/")
mycounts<-read.table("Mergedata.txt",header = TRUE,row.names = 1,sep ="\t" )
mycounts_new<-mycounts
a<-colnames(mycounts)
c<-data.frame(0,0,0,0,0,0,0,0)
colnames(c)<-c("j","k","l","p","q","r","y","z")

for (j in 1:6) {
  for (k in (j+1):7) {
    
    for (l in (k+1):8) {
      
      for (p in 9:17) {
        for (q in (p+1):18) {
          
          for (r in (q+1):19) {
            
            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r])]
            condition<-factor(c(rep("KO",3),rep("WT",3)),levels = c("WT","KO"))
            colData<-data.frame(row.names = colnames(mycounts),condition)
            
            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
            dds <- DESeq(dds)
            
            res= as.data.frame(results(dds))
            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
            y<-length(rownames(res_hi))
            z<-length(rownames(res_lo))
            if(y>z & (y/z)>=4){
              
              b<-as.data.frame(t(as.data.frame(c(j,k,l,p,q,r,y,z))))
              colnames(b)<-c("j","k","l","p","q","r","y","z")
              c<-rbind(c,b)
              write.table(mycounts,file = paste0("mycounts_new",j,k,l,p,q,r,".txt"),sep = "\t")
            }
            
          }
        }
      }
    }
  }
}