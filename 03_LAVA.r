library(LAVA)
library(dplyr)

input = process.input(input.info.file="/storage/zhenghoufengLab/qianyu/project/RA_CVD/input.txt",           # input info file
                      sample.overlap.file="/storage/zhenghoufengLab/qianyu/project/RA_CVD/sample.overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="/storage/zhenghoufengLab/qianyu/Software/g1000_eur",                    # reference genotype data prefix
                      phenos=c("RA","AF"))       # subset of phenotypes listed in the input info file that we want to process

loci = read.loci("/storage/zhenghoufengLab/qianyu/Software/LAVA/Locus.txt")

newdata_1=matrix(,ncol=4)
colnames(newdata_1) <- c("phen","h2.obs","p","loci")

newdata_2=matrix(,ncol=10)
colnames(newdata_2) <- c("phen1","phen2","rho","rho.lower","rho.upper","r2","r2.lower","r2.upper","p","loci")

nrow(loci)

for(n in 1:nrow(loci)){
locus = process.locus(loci[n,], input)
if(!is.null(locus)) {
    test<-run.univ(locus)
    if(nrow(test)>1){
        test_2<-filter(.data=test,phen=="RA"&p<0.05)
        if(nrow(test_2)==1){
            test_3<-filter(.data=test,p<0.05)              
                if(nrow(test_3)>1){
loc.out_2=run.univ.bivar(locus,target="RA")
tem_1<-loc.out_2$univ
tem_1$loci<-n
newdata_1<-rbind(newdata_1,tem_1)

tem_2<-loc.out_2$bivar
tem_2$loci<-n
newdata_2<-rbind(newdata_2,tem_2)
                }
        }
}
}
print(n)
}

write.table(newdata_1,"/storage/zhenghoufengLab/qianyu/project/RA_CVD/LAVA/RA_AF_LAVA_univ_samplecor.txt",sep="\t",quote=F)
write.table(newdata_2,"/storage/zhenghoufengLab/qianyu/project/RA_CVD/LAVA/RA_AF_LAVA_bivar_samplecor.txt",sep="\t",quote=F)
