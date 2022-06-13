library(tidyverse)
library(DESeq2)
#import data
setwd("D:/yyy/")
mycounts<-read.table("Mergedata.txt",header = TRUE,row.names = 1,sep ="\t" )
mycounts_new<-mycounts
a<-colnames(mycounts)
c<-data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
colnames(c)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
for (i in 2:5) {
  for (j in 1:1) {
    for (k in (j+1):2) {
      if(i>2 ){
        for (l in (k+1):3) {
          if(i>3 ){
            for (m in (l+1):4) {
              if(i>4 ){
                for (n in (m+1):5) {
                  
                  for (o in 2:9) {
                    for (p in 6:13) {
                      for (q in (p+1):7) {
                        if(o>2 ){
                          for (r in (q+1):8) {
                            if(o>3 ){
                              for (s in (r+1):9) {
                                if(o>4 ){
                                  for (t in (s+1):10) {
                                    if(o>5 ){
                                      for (u in (t+1):11) {
                                        if(o>6 ){
                                          for (v in (u+1):12) {
                                            if(o>7 ){
                                              for (w in (v+1):13) {
                                                if(o>8 ){
                                                  for (x in (w+1):14) {
                                                    
                                                    mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w],a[x])]
                                                    condition<-factor(c(rep("KO",5),rep("WT",9)),levels = c("WT","KO"))
                                                    colData<-data.frame(row.names = colnames(mycounts),condition)
                                                    
                                                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                                    dds <- DESeq(dds)
                                                    
                                                    res= as.data.frame(results(dds))
                                                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                                    y<-length(rownames(res_hi))
                                                    z<-length(rownames(res_lo))
                                                    if(y>z){
                                                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                                      c<-rbind(c,b)
                                                    }
                                                  }
                                                  }
                                                }
                                                mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w])]
                                                condition<-factor(c(rep("KO",5),rep("WT",8)),levels = c("WT","KO"))
                                                colData<-data.frame(row.names = colnames(mycounts),condition)
                                                
                                                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                                dds <- DESeq(dds)
                                                
                                                res= as.data.frame(results(dds))
                                                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                                y<-length(rownames(res_hi))
                                                z<-length(rownames(res_lo))
                                                if(y>z){
                                                  x<-0
                                                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                                  c<-rbind(c,b)
                                                }
                                                
                                                }
                                            }
                                            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s],a[t],a[u],a[v])]
                                            condition<-factor(c(rep("KO",5),rep("WT",7)),levels = c("WT","KO"))
                                            colData<-data.frame(row.names = colnames(mycounts),condition)
                                            
                                            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                            dds <- DESeq(dds)
                                            
                                            res= as.data.frame(results(dds))
                                            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                            y<-length(rownames(res_hi))
                                            z<-length(rownames(res_lo))
                                            if(y>z){
                                              w<-0
                                              x<-0
                                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                              c<-rbind(c,b)
                                            }
                                            
                                            }
                                        }
                                        mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s],a[t],a[u])]
                                        condition<-factor(c(rep("KO",5),rep("WT",6)),levels = c("WT","KO"))
                                        colData<-data.frame(row.names = colnames(mycounts),condition)
                                        
                                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                        dds <- DESeq(dds)
                                        
                                        res= as.data.frame(results(dds))
                                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                        y<-length(rownames(res_hi))
                                        z<-length(rownames(res_lo))
                                        if(y>z){
                                          v<-0
                                          w<-0
                                          x<-0
                                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                          c<-rbind(c,b)
                                        }
                                        
                                        }
                                    }
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
                                      u<-0
                                      v<-0
                                      w<-0
                                      x<-0
                                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                      c<-rbind(c,b)
                                    }
                                    
                                    }
                                }
                                mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s])]
                                condition<-factor(c(rep("KO",5),rep("WT",4)),levels = c("WT","KO"))
                                colData<-data.frame(row.names = colnames(mycounts),condition)
                                
                                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                dds <- DESeq(dds)
                                
                                res= as.data.frame(results(dds))
                                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                y<-length(rownames(res_hi))
                                z<-length(rownames(res_lo))
                                if(y>z){
                                  t<-0
                                  u<-0
                                  v<-0
                                  w<-0
                                  x<-0
                                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                  c<-rbind(c,b)
                                }
                                
                                }
                            }
                            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r])]
                            condition<-factor(c(rep("KO",5),rep("WT",3)),levels = c("WT","KO"))
                            colData<-data.frame(row.names = colnames(mycounts),condition)
                            
                            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                            dds <- DESeq(dds)
                            
                            res= as.data.frame(results(dds))
                            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                            y<-length(rownames(res_hi))
                            z<-length(rownames(res_lo))
                            if(y>z){
                              s<-0
                              t<-0
                              u<-0
                              v<-0
                              w<-0
                              x<-0
                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                              c<-rbind(c,b)
                            }
                            
                            }
                        }
                        mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q])]
                        condition<-factor(c(rep("KO",5),rep("WT",2)),levels = c("WT","KO"))
                        colData<-data.frame(row.names = colnames(mycounts),condition)
                        
                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                        dds <- DESeq(dds)
                        
                        res= as.data.frame(results(dds))
                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                        y<-length(rownames(res_hi))
                        z<-length(rownames(res_lo))
                        if(y>z){
                          r<-0
                          s<-0
                          t<-0
                          u<-0
                          v<-0
                          w<-0
                          x<-0
                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                          c<-rbind(c,b)
                        }
                        
                        }
                    }
                  }
                }
              }
              for (o in 2:9) {
                for (p in 6:13) {
                  for (q in (p+1):14) {
                    if(o>2 & p<13){
                      for (r in (q+1):14) {
                        if(o>3 & q<13){
                          for (s in (r+1):14) {
                            if(o>4 & r <13){
                              for (t in (s+1):14) {
                                if(o>5 & s<13){
                                  for (u in (t+1):14) {
                                    if(o>6 & t<13){
                                      for (v in (u+1):14) {
                                        if(o>7 & u<13){
                                          for (w in (v+1):14) {
                                            if(o>8 & v<13){
                                              for (x in (w+1):14) {
                                                if(w<13){
                                                mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w],a[x])]
                                                condition<-factor(c(rep("KO",4),rep("WT",9)),levels = c("WT","KO"))
                                                colData<-data.frame(row.names = colnames(mycounts),condition)
                                                
                                                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                                dds <- DESeq(dds)
                                                
                                                res= as.data.frame(results(dds))
                                                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                                y<-length(rownames(res_hi))
                                                z<-length(rownames(res_lo))
                                                if(y>z){
                                                  n<-0
                                                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                                  c<-rbind(c,b)
                                                }
                                              }
                                              }
                                            }
                                            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w])]
                                            condition<-factor(c(rep("KO",4),rep("WT",8)),levels = c("WT","KO"))
                                            colData<-data.frame(row.names = colnames(mycounts),condition)
                                            
                                            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                            dds <- DESeq(dds)
                                            
                                            res= as.data.frame(results(dds))
                                            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                            y<-length(rownames(res_hi))
                                            z<-length(rownames(res_lo))
                                            if(y>z){
                                              n<-0
                                              x<-0
                                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                              c<-rbind(c,b)
                                            }
                                            
                                          }
                                        }
                                        mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r],a[s],a[t],a[u],a[v])]
                                        condition<-factor(c(rep("KO",4),rep("WT",7)),levels = c("WT","KO"))
                                        colData<-data.frame(row.names = colnames(mycounts),condition)
                                        
                                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                        dds <- DESeq(dds)
                                        
                                        res= as.data.frame(results(dds))
                                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                        y<-length(rownames(res_hi))
                                        z<-length(rownames(res_lo))
                                        if(y>z){
                                          n<-0
                                          w<-0
                                          x<-0
                                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                          c<-rbind(c,b)
                                        }
                                        
                                      }
                                    }
                                    mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r],a[s],a[t],a[u])]
                                    condition<-factor(c(rep("KO",4),rep("WT",6)),levels = c("WT","KO"))
                                    colData<-data.frame(row.names = colnames(mycounts),condition)
                                    
                                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                    dds <- DESeq(dds)
                                    
                                    res= as.data.frame(results(dds))
                                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                    y<-length(rownames(res_hi))
                                    z<-length(rownames(res_lo))
                                    if(y>z){
                                      n<-0
                                      v<-0
                                      w<-0
                                      x<-0
                                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                      c<-rbind(c,b)
                                    }
                                    
                                  }
                                }
                                mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r],a[s],a[t])]
                                condition<-factor(c(rep("KO",4),rep("WT",5)),levels = c("WT","KO"))
                                colData<-data.frame(row.names = colnames(mycounts),condition)
                                
                                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                dds <- DESeq(dds)
                                
                                res= as.data.frame(results(dds))
                                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                y<-length(rownames(res_hi))
                                z<-length(rownames(res_lo))
                                if(y>z){
                                  n<-0
                                  u<-0
                                  v<-0
                                  w<-0
                                  x<-0
                                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                  c<-rbind(c,b)
                                }
                                
                              }
                            }
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
                            if(y>z){
                              n<-0
                              t<-0
                              u<-0
                              v<-0
                              w<-0
                              x<-0
                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                              c<-rbind(c,b)
                            }
                            
                          }
                        }
                        mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q],a[r])]
                        condition<-factor(c(rep("KO",4),rep("WT",3)),levels = c("WT","KO"))
                        colData<-data.frame(row.names = colnames(mycounts),condition)
                        
                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                        dds <- DESeq(dds)
                        
                        res= as.data.frame(results(dds))
                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                        y<-length(rownames(res_hi))
                        z<-length(rownames(res_lo))
                        if(y>z){
                          n<-0
                          s<-0
                          t<-0
                          u<-0
                          v<-0
                          w<-0
                          x<-0
                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                          c<-rbind(c,b)
                        }
                        
                      }
                    }
                    mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[p],a[q])]
                    condition<-factor(c(rep("KO",4),rep("WT",2)),levels = c("WT","KO"))
                    colData<-data.frame(row.names = colnames(mycounts),condition)
                    
                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                    dds <- DESeq(dds)
                    
                    res= as.data.frame(results(dds))
                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                    y<-length(rownames(res_hi))
                    z<-length(rownames(res_lo))
                    if(y>z){
                      n<-0
                      r<-0
                      s<-0
                      t<-0
                      u<-0
                      v<-0
                      w<-0
                      x<-0
                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                      c<-rbind(c,b)
                    }
                    
                  }
                }
              }
              }
          }
          for (o in 2:9) {
            for (p in 6:13) {
              for (q in (p+1):14) {
                if(o>2 & p<13){
                  for (r in (q+1):14) {
                    if(o>3 & q<13){
                      for (s in (r+1):14) {
                        if(o>4 & r <13){
                          for (t in (s+1):14) {
                            if(o>5 & s<13){
                              for (u in (t+1):14) {
                                if(o>6 & t<13){
                                  for (v in (u+1):14) {
                                    if(o>7 & u<13){
                                      for (w in (v+1):14) {
                                        if(o>8 & v<13){
                                          for (x in (w+1):14) {
                                            if(w<13){
                                            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w],a[x])]
                                            condition<-factor(c(rep("KO",3),rep("WT",9)),levels = c("WT","KO"))
                                            colData<-data.frame(row.names = colnames(mycounts),condition)
                                            
                                            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                            dds <- DESeq(dds)
                                            
                                            res= as.data.frame(results(dds))
                                            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                            y<-length(rownames(res_hi))
                                            z<-length(rownames(res_lo))
                                            if(y>z){
                                              m<-0
                                              n<-0
                                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                              c<-rbind(c,b)
                                            }
                                          }
                                          }
                                        }
                                        mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w])]
                                        condition<-factor(c(rep("KO",3),rep("WT",8)),levels = c("WT","KO"))
                                        colData<-data.frame(row.names = colnames(mycounts),condition)
                                        
                                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                        dds <- DESeq(dds)
                                        
                                        res= as.data.frame(results(dds))
                                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                        y<-length(rownames(res_hi))
                                        z<-length(rownames(res_lo))
                                        if(y>z){
                                          m<-0
                                          n<-0
                                          x<-0
                                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                          c<-rbind(c,b)
                                        }
                                        
                                      }
                                    }
                                    mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r],a[s],a[t],a[u],a[v])]
                                    condition<-factor(c(rep("KO",3),rep("WT",7)),levels = c("WT","KO"))
                                    colData<-data.frame(row.names = colnames(mycounts),condition)
                                    
                                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                    dds <- DESeq(dds)
                                    
                                    res= as.data.frame(results(dds))
                                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                    y<-length(rownames(res_hi))
                                    z<-length(rownames(res_lo))
                                    if(y>z){
                                      m<-0
                                      n<-0
                                      w<-0
                                      x<-0
                                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                      c<-rbind(c,b)
                                    }
                                    
                                  }
                                }
                                mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r],a[s],a[t],a[u])]
                                condition<-factor(c(rep("KO",3),rep("WT",6)),levels = c("WT","KO"))
                                colData<-data.frame(row.names = colnames(mycounts),condition)
                                
                                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                dds <- DESeq(dds)
                                
                                res= as.data.frame(results(dds))
                                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                y<-length(rownames(res_hi))
                                z<-length(rownames(res_lo))
                                if(y>z){
                                  m<-0
                                  n<-0
                                  v<-0
                                  w<-0
                                  x<-0
                                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                  c<-rbind(c,b)
                                }
                                
                              }
                            }
                            mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r],a[s],a[t])]
                            condition<-factor(c(rep("KO",3),rep("WT",5)),levels = c("WT","KO"))
                            colData<-data.frame(row.names = colnames(mycounts),condition)
                            
                            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                            dds <- DESeq(dds)
                            
                            res= as.data.frame(results(dds))
                            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                            y<-length(rownames(res_hi))
                            z<-length(rownames(res_lo))
                            if(y>z){
                              m<-0
                              n<-0
                              u<-0
                              v<-0
                              w<-0
                              x<-0
                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                              c<-rbind(c,b)
                            }
                            
                          }
                        }
                        mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q],a[r],a[s])]
                        condition<-factor(c(rep("KO",3),rep("WT",4)),levels = c("WT","KO"))
                        colData<-data.frame(row.names = colnames(mycounts),condition)
                        
                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                        dds <- DESeq(dds)
                        
                        res= as.data.frame(results(dds))
                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                        y<-length(rownames(res_hi))
                        z<-length(rownames(res_lo))
                        if(y>z){
                          m<-0
                          n<-0
                          t<-0
                          u<-0
                          v<-0
                          w<-0
                          x<-0
                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                          c<-rbind(c,b)
                        }
                        
                      }
                    }
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
                    if(y>z){
                      m<-0
                      n<-0
                      s<-0
                      t<-0
                      u<-0
                      v<-0
                      w<-0
                      x<-0
                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                      c<-rbind(c,b)
                    }
                    
                  }
                }
                mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[p],a[q])]
                condition<-factor(c(rep("KO",3),rep("WT",2)),levels = c("WT","KO"))
                colData<-data.frame(row.names = colnames(mycounts),condition)
                
                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                dds <- DESeq(dds)
                
                res= as.data.frame(results(dds))
                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                y<-length(rownames(res_hi))
                z<-length(rownames(res_lo))
                if(y>z){
                  m<-0
                  n<-0
                  r<-0
                  s<-0
                  t<-0
                  u<-0
                  v<-0
                  w<-0
                  x<-0
                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                  c<-rbind(c,b)
                }
                
              }
            }
          }
          }
      }
      for (o in 2:9) {
        for (p in 6:13) {
          for (q in (p+1):14) {
            if(o>2 & p<13){
              for (r in (q+1):14) {
                if(o>3 & q<13){
                  for (s in (r+1):14) {
                    if(o>4 & r <13){
                      for (t in (s+1):14) {
                        if(o>5 & s<13){
                          for (u in (t+1):14) {
                            if(o>6 & t<13){
                              for (v in (u+1):14) {
                                if(o>7 & u<13){
                                  for (w in (v+1):14) {
                                    if(o>8 & v<13){
                                      for (x in (w+1):14) {
                                        if(w<13){
                                        mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w],a[x])]
                                        condition<-factor(c(rep("KO",2),rep("WT",9)),levels = c("WT","KO"))
                                        colData<-data.frame(row.names = colnames(mycounts),condition)
                                        
                                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                        dds <- DESeq(dds)
                                        
                                        res= as.data.frame(results(dds))
                                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                        y<-length(rownames(res_hi))
                                        z<-length(rownames(res_lo))
                                        if(y>z){
                                          l<-0
                                          m<-0
                                          n<-0
                                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                          c<-rbind(c,b)
                                        }
                                      }
                                      }
                                    }
                                    mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q],a[r],a[s],a[t],a[u],a[v],a[w])]
                                    condition<-factor(c(rep("KO",2),rep("WT",8)),levels = c("WT","KO"))
                                    colData<-data.frame(row.names = colnames(mycounts),condition)
                                    
                                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                    dds <- DESeq(dds)
                                    
                                    res= as.data.frame(results(dds))
                                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                    y<-length(rownames(res_hi))
                                    z<-length(rownames(res_lo))
                                    if(y>z){
                                      l<-0
                                      m<-0
                                      n<-0
                                      x<-0
                                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                      c<-rbind(c,b)
                                    }
                                    
                                  }
                                }
                                mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q],a[r],a[s],a[t],a[u],a[v])]
                                condition<-factor(c(rep("KO",2),rep("WT",7)),levels = c("WT","KO"))
                                colData<-data.frame(row.names = colnames(mycounts),condition)
                                
                                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                                dds <- DESeq(dds)
                                
                                res= as.data.frame(results(dds))
                                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                                y<-length(rownames(res_hi))
                                z<-length(rownames(res_lo))
                                if(y>z){
                                  l<-0
                                  m<-0
                                  n<-0
                                  w<-0
                                  x<-0
                                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                                  c<-rbind(c,b)
                                }
                                
                              }
                            }
                            mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q],a[r],a[s],a[t],a[u])]
                            condition<-factor(c(rep("KO",2),rep("WT",6)),levels = c("WT","KO"))
                            colData<-data.frame(row.names = colnames(mycounts),condition)
                            
                            dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                            dds <- DESeq(dds)
                            
                            res= as.data.frame(results(dds))
                            res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                            res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                            y<-length(rownames(res_hi))
                            z<-length(rownames(res_lo))
                            if(y>z){
                              l<-0
                              m<-0
                              n<-0
                              v<-0
                              w<-0
                              x<-0
                              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                              c<-rbind(c,b)
                            }
                            
                          }
                        }
                        mycounts<-mycounts_new[,c(a[j],a[k],a[n],a[p],a[q],a[r],a[s],a[t])]
                        condition<-factor(c(rep("KO",2),rep("WT",5)),levels = c("WT","KO"))
                        colData<-data.frame(row.names = colnames(mycounts),condition)
                        
                        dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                        dds <- DESeq(dds)
                        
                        res= as.data.frame(results(dds))
                        res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                        res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                        y<-length(rownames(res_hi))
                        z<-length(rownames(res_lo))
                        if(y>z){
                          l<-0
                          m<-0
                          n<-0
                          u<-0
                          v<-0
                          w<-0
                          x<-0
                          b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                          colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                          c<-rbind(c,b)
                        }
                        
                      }
                    }
                    mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q],a[r],a[s])]
                    condition<-factor(c(rep("KO",2),rep("WT",4)),levels = c("WT","KO"))
                    colData<-data.frame(row.names = colnames(mycounts),condition)
                    
                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                    dds <- DESeq(dds)
                    
                    res= as.data.frame(results(dds))
                    res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                    res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                    y<-length(rownames(res_hi))
                    z<-length(rownames(res_lo))
                    if(y>z){
                      l<-0
                      m<-0
                      n<-0
                      t<-0
                      u<-0
                      v<-0
                      w<-0
                      x<-0
                      b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                      colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                      c<-rbind(c,b)
                    }
                    
                  }
                }
                mycounts<-mycounts_new[,c(a[j],a[k],a[p],a[q],a[r])]
                condition<-factor(c(rep("KO",2),rep("WT",3)),levels = c("WT","KO"))
                colData<-data.frame(row.names = colnames(mycounts),condition)
                
                dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                dds <- DESeq(dds)
                
                res= as.data.frame(results(dds))
                res_hi<-res[which(res$padj <0.05 & res$log2FoldChange >0.58),]
                res_lo<-res[which(res$padj <0.05 & res$log2FoldChange <(-0.58)),]
                y<-length(rownames(res_hi))
                z<-length(rownames(res_lo))
                if(y>z){
                  l<-0
                  m<-0
                  n<-0
                  s<-0
                  t<-0
                  u<-0
                  v<-0
                  w<-0
                  x<-0
                  b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
                  colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
                  c<-rbind(c,b)
                }
                
              }
            }
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
            if(y>z){
              l<-0
              m<-0
              n<-0
              r<-0
              s<-0
              t<-0
              u<-0
              v<-0
              w<-0
              x<-0
              b<-as.data.frame(t(as.data.frame(c(j,k,l,m,n,p,q,r,s,t,u,v,w,x,y,z))))
              colnames(b)<-c("j","k","l","m","n","p","q","r","s","t","u","v","w","x","y","z")
              c<-rbind(c,b)
            }
            
          }
        }
      }
      }
  }
}
}



