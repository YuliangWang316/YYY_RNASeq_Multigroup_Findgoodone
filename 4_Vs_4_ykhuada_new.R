library(tidyverse)
library(DESeq2)
#import data
setwd("D:/yyy/Merge_HUADA_YK_HUADA2/")
mycounts<-read.table("Mergedata_FOB.txt",header = TRUE,row.names = 1,sep ="\t" )
mycounts_new<-mycounts
a<-colnames(mycounts)
c<-data.frame(0,0,0,0,0,0,0,0,0,0)
colnames(c)<-c("j","k","l","m","p","q","r","s","y","z")

for (j in 1:7) {
  for (k in (j+1):8) {
    
    for (l in (k+1):9) {
      for (m in (l+1):10) {
        
      
      for (p in 11:24) {
        for (q in (p+1):25) {
          
          for (r in (q+1):26) {
            for (s in (r+1):27) {
              
            
            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r],a[s])]
            condition<-factor(c(rep("KO",4),rep("WT",4)),levels = c("WT","KO"))
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
            if(y>z & (y/z) >= 4 & prdm1_b<0.05 & prdm1_a>0.58 & is.na(prdm1_b) == FALSE & y>=100){
              
              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,p,q,r,s,y,z))))
              colnames(b)<-c("j","k","l","m","p","q","r","s","y","z")
              c<-rbind(c,b)
              write.table(mycounts,file = paste0("mycounts_YKHUADA_4vs4_new",j,k,l,m,p,q,r,s,y,z,prdm1_a,".txt"),sep = "\t")
              write.table(res,file = paste0("mycounts_YKHUADA_Diff_4vs4_new",j,k,l,m,p,q,r,s,y,z,prdm1_a,".txt"),sep = "\t")
            }
            
          }
        }
      }
    }
  }
}
  }
}
