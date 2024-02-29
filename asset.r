library("ASSET")
library(data.table)

CVD<-c("AF","CAD","AS","HF")

N00N11N01<-read.table("/storage/zhenghoufengLab/qianyu/project/RA_CVD/ASSET/N00N11N10.txt",head=T)
N00N11N01$RA<-as.numeric(N00N11N01$RA)
N00N11N01$AF<-as.numeric(N00N11N01$AF)

for(n in 1:1){
    m=2*(n-1)+1
    data<-read.table(paste("/storage/zhenghoufengLab/qianyu/project/RA_CVD/ASSET/RA_",CVD[n],"_meta.csv",sep=""),sep=",",head=TRUE)
    data<-data[,c(1,2,4,3,5)]
    snps<-as.vector(data[,"SNP"])
    traits.lab<-c("RA",CVD[n])
    beta.hat<-as.matrix(data[, paste("beta_",traits.lab,sep="")])
    sigma.hat<-as.matrix(data[, paste("se_",traits.lab,sep="")])
    
    N00<-as.matrix(N00N11N01[c(m,m+1),c(2,3)])
    colnames(N00)<-c("RA",CVD[n])
    rownames(N00)<-c("RA",CVD[n])
    
    N11<-as.matrix(N00N11N01[c(m+8,m+9),c(2,3)])
    colnames(N11)<-c("RA",CVD[n])
    rownames(N11)<-c("RA",CVD[n]) 
    
    N10<-as.matrix(N00N11N01[c(m+16,m+17),c(2,3)])
    colnames(N10)<-c("RA",CVD[n])
    rownames(N10)<-c("RA",CVD[n])

    cor<-list(N11=N11, N00=N00, N10=N10)

    ncase<-diag(N11)
    ncntl<-diag(N00)
    #ncase<-c(N00N11N01[c(n+8),2],N00N11N01[c(n+9),3])
    #ncntl<-c(N00N11N01[c(n),2],N00N11N01[c(n+1),3])

 # Now let us call h.traits on these summary data. 
 res <- h.traits(snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=cor, 
                 cor.numr=FALSE, search=NULL, side=2, meta=TRUE, 
                 zmax.args=NULL, meth.pval="DLM")
newdata<-cbind(res[["Subset.2sided"]][["beta.1"]],res[["Subset.2sided"]][["sd.1"]],res[["Subset.2sided"]][["pval.1"]],res[["Subset.2sided"]][["beta.2"]],res[["Subset.2sided"]][["sd.2"]],res[["Subset.2sided"]][["pval.2"]],res[["Subset.2sided"]][["pval"]],res[["Subset.2sided"]][["pheno.1"]],res[["Subset.2sided"]][["pheno.2"]])
colnames(newdata)<-c("beta_1","sd_1","P_1","beta_2","sd_2","P_2","P_meta","Trait_1_pos","Trait_2_pos","Trait_1_neg","Trait_2_neg")
write.table(newdata,paste("/storage/zhenghoufengLab/qianyu/project/RA_CVD/ASSET/ASSET_RA",CVD[n],".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
}
