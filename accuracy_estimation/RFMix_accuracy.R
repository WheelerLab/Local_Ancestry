###calc RFMix accuracy

library(data.table)
library(argparse)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

parser <- ArgumentParser()
parser$add_argument("--viterbi", help="Rdata object ouput by MOSAIC")
parser$add_argument("--haps.hap.gz", help="admixed sample list")
parser$add_argument("--haps.sample", help="admixed samle list")
parser$add_argument("--nanc", help="number of ancestries estimated")
parser$add_argument("--result", help="results file output by adsim")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

print("processing snp ids")
snps<-fread("zcat " %&% args$haps.hap.gz, select = c(1,3))
colnames(snps)<-c("chm","pos")
rfout<-fread(args$viterbi, header = F)
rfout<-as.data.frame(cbind.data.frame(snps,rfout))
true_ancestry<-fread(args$result, header = T)
true_ancestry_subset<-inner_join(true_ancestry,snps,by=c("chm","pos"))

# str(true_ancestry_subset)
# str(rfout)
##calculate haploid correlation
hap_cor<-c(rep(NA,40))
for (i in c(1:40)){
  x<-select(true_ancestry_subset, (i+2)) %>% unlist %>% as.numeric()
  y<-select(rfout,(i+2)) %>% unlist %>% as.numeric()
  cor<-cor.test(x,y,method="pearson")
  hap_cor[i]<-cor$estimate
}
hap_cor[is.na(hap_cor)]<-1.01
hap_cor<-as.list(hap_cor)
fwrite(hap_cor,args$out %&% "_haps_accuracy.txt", sep="\t", col.names = F)
#instantiate empty df
#convert haplod ancestry to diploid
diploid_true<-data.frame(matrix(ncol=20,nrow=dim(true_ancestry_subset)[1]))
diploid_est<-data.frame(matrix(ncol=20,nrow=dim(rfout)[1]))

if (args$nanc == 3){
  random<-c(1,2,3)
  for ( i in c(1:20)){ 
    cat( i, "/20\n")
    hap2<-i*2
    hap1<-hap2-1
    dip_est<-paste(rfout[,(hap1+2)],rfout[,(hap2+2)],sep="")
    dip_est<-gsub("32","23",gsub("31","13",gsub("21","12",dip_est)))
    dip_est[dip_est == "11"] <- 1
    dip_est[dip_est == "12"] <- 2
    dip_est[dip_est == "13"] <- 3
    dip_est[dip_est == "22"] <- 4
    dip_est[dip_est == "23"] <- 5
    dip_est[dip_est == "33"] <- 6
    dip_est<-as.numeric(dip_est)
    diploid_est[,i]<-dip_est
    
    #repeat process for results df
    dip_t<-paste(true_ancestry_subset[,hap1+2],true_ancestry_subset[,hap2+2],sep="")
    dip_t<-gsub("32","23",gsub("31","13",gsub("21","12",dip_t)))
    dip_t[dip_t == "11"] <- 1
    dip_t[dip_t == "12"] <- 2
    dip_t[dip_t == "13"] <- 3
    dip_t[dip_t == "22"] <- 4
    dip_t[dip_t == "23"] <- 5
    dip_t[dip_t == "33"] <- 6
    dip_t<-as.numeric(dip_t)
    diploid_true[,i]<-dip_t
  }
} else if (args$nanc ==2){
  random<-c(1,2)
  for ( i in c(1:20)){ 
    cat( i, "/20\n")
    hap2<-i*2 
    hap1<-hap2-1
    dip_est<-paste(rfout[,(hap1+2)],rfout[,(hap2+2)],sep="")
    dip_est<-gsub("21","12",dip_est)
    dip_est[dip_est == "11"] <- 1
    dip_est[dip_est == "12"] <- 2
    dip_est[dip_est == "22"] <- 3
    dip_est<-as.numeric(dip_est)
    diploid_est[,i]<-dip_est
    
    #repeat process for results df
    dip_t<-paste(true_ancestry_subset[,hap1+2],true_ancestry_subset[,hap2+2],sep="")
    dip_t<-gsub("21","12",dip_t)
    dip_t[dip_t == "11"] <- 1
    dip_t[dip_t == "12"] <- 2
    dip_t[dip_t == "22"] <- 3
    dip_t<-as.numeric(dip_t)
    diploid_true[,i]<-dip_t
  }
}

dip_cor<-c(rep(NA,20))
for (i in c(1:20)){
  x<-diploid_true[,i]
  y<-diploid_est[,i]
  dcor<-cor.test(x,y,method="pearson")
  dip_cor[i]<-dcor$estimate
}
dip_cor[is.na(dip_cor)]<-1.01
dip_cor<-as.list(dip_cor)
fwrite(dip_cor,args$out %&% "_dip_accuracy.txt",col.names = F, sep ='\t')