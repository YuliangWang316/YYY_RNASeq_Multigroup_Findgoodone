library(tidyverse)
library(DESeq2)
#import data
setwd("D:/yyy/")
mycounts<-read.table("Mergedata.txt",header = TRUE,row.names = 1,sep ="\t" )
mycounts_new<-mycounts
a<-colnames(mycounts)
c<-data.frame(0,0,0,0,0,0,0,0,0,0,0,0)
colnames(c)<-c("j","k","l","m","n","p","q","r","s","t","y","z")

  for (j in 1:1) {
    for (k in (j+1):2) {
      
        for (l in (k+1):3) {
          
            for (m in (l+1):4) {
              
                for (n in (m+1):5) {
                  
                    for (p in 6:10) {
                      for (q in (p+1):11) {
                       
                          for (r in (q+1):12) {
                            
                              for (s in (r+1):13) {
                                
                                  for (t in (s+1):14) {
                                    
                                    mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s],a[t])]
                                    condition<-factor(c(rep("KO",5),rep("WT",5)),levels = c("WT","KO"))
                                    colData<-data.frame(row.names = colnames(mycounts),condition)
                                    
                                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                    dds <- DESeq(dds)
                                    
                                    res= as.data.frame(results(dds))
                                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                    y<-length(rownames(res_hi))
                                    z<-length(rownames(res_lo))
                                    if(y>z){
                                      
                                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,y,z))))
                                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","y","z")
                                      c<-rbind(c,b)
                                    }
                                    
                                    }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                }
              }




