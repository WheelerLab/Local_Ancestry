###calc LAMP accuracy

library(data.table)
library(argparse)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

parser <- ArgumentParser()
parser$add_argument("--long", help="Rdata object ouput by MOSAIC")
parser$add_argument("--haps.hap.gz", help="admixed sample list")
parser$add_argument("--haps.sample", help="admixed samle list")
parser$add_argument("--result", help="results file output by adsim")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

print("processing sample ids")
##read in samples
samps<-fread(args$haps.sample, header = F, skip = 2)
samps$V1<-paste(samps$V1, ".0", sep = "")
samps$V2<-paste(samps$V2, ".1", sep = "")
samps$V3<-NULL
#assign sample ids to MOSAIC
ids<-as.vector(t(samps))

print("processing snp ids")
##read in snps
snps<-fread("zcat " %&% args$haps.hap.gz, select = c(1,3))
colnames(snps)<-c("chm","pos")

print("processing true ancestry")
res<-fread(args$result, header = T) %>% inner_join(snps, by = c("chm","pos")) 

print("processing in LAMPLD results")
long<-fread(args$long, header = F)
long<-strsplit(long$V1,"")
long<-matrix(unlist(long), nrow=length(long), byrow=T)  %>% t()
long<-apply(long,2,as.numeric) %>% as.data.frame(stringsAsFactors = F) 
diploid_true<-data.frame(matrix(ncol = 20,nrow=50000))
diploid_predicted<-data.frame(matrix(ncol = 20,nrow=50000))
# str(long)
print("calculating haploid accuracy")
## calc both true and maximized
hap_cor_unmaxed<-c(rep(NA,40))
hap_cor_maximized<-c(rep(NA,40))
for (i in c(1:20)){
  hap2<- i * 2
  hap1<-hap2-1
  x1<-select(long, (hap1)) %>% unlist %>% as.numeric()
  y1<-select(res,(hap1)) %>% unlist %>% as.numeric()
  
  x2<-select(long, (hap2)) %>% unlist %>% as.numeric()
  y2<-select(res,(hap2)) %>% unlist %>% as.numeric()
  
  hap1_unflipped<-cor.test(x1,y1,method="pearson")
  hap2_unflipped<-cor.test(x2,y2,method="pearson")
  hap_cor_unmaxed[hap1]<-hap1_unflipped$estimate
  hap_cor_unmaxed[hap2]<-hap2_unflipped$estimate
  
  hap1_flipped<-cor.test(x1,y2,method="pearson")
  hap2_flipped<-cor.test(x2,y1,method="pearson")
  # str(hap1_unflipped$estimate)
  # str(hap2_unflipped$estimate)
  # str(hap1_flipped$estimate)
  # str(hap2_flipped$estimate)
  mu_unflipped<-mean(c(hap1_unflipped$estimate,hap2_unflipped$estimate),na.rm = T)
  mu_flipped<-mean(c(hap1_flipped$estimate,hap2_flipped$estimate), na.rm = T)
  mu_unflipped<-ifelse(!is.na(mu_unflipped),mu_unflipped,0) 
  mu_flipped<-ifelse(!is.na(mu_flipped),mu_flipped,0) 
  # print("here")
  if (mu_unflipped > mu_flipped){
    hap_cor_maximized[hap1]<-ifelse(!is.na(hap1_unflipped$estimate),hap1_unflipped$estimate,1.01)
    hap_cor_maximized[hap2]<-ifelse(!is.na(hap2_unflipped$estimate),hap2_unflipped$estimate,1.01)
  } else {
    hap_cor_maximized[hap1]<-ifelse(!is.na(hap1_flipped$estimate),hap1_unflipped$estimate,1.01)
    hap_cor_maximized[hap2]<-ifelse(!is.na(hap2_flipped$estimate),hap2_unflipped$estimate,1.01)
  }
    
}
hap_cor_unmaxed[is.na(hap_cor_unmaxed)]<-1.01
hap_cor_maximized[is.na(hap_cor_maximized)]<-1.01
hap_cor_unmaxed<-as.list(hap_cor_unmaxed)
hap_cor_maximized<-as.list(hap_cor_maximized)
str(hap_cor_unmaxed)
str(hap_cor_maximized)
fwrite(hap_cor_unmaxed,args$out %&% "_haps_unmaxed_accuracy.txt", sep="\t", col.names = F)
fwrite(hap_cor_maximized,args$out %&% "_haps_maxed_accuracy.txt", sep="\t", col.names = F)

print("coverting ancestries to dosage")
for (i in c(1:20)){
  cat(i,"/ 20\n")
  hap2<- i * 2
  hap1<-hap2-1
  
  dip_pred<-paste(long[,hap1],long[,hap2],sep="")
  dip_pred[dip_pred == "00"] <- 1
  dip_pred[dip_pred == "01"] <- 2
  dip_pred[dip_pred == "02"] <- 3
  dip_pred[dip_pred == "11"] <- 4
  dip_pred[dip_pred == "12"] <- 5
  dip_pred[dip_pred == "22"] <- 6
  dip_pred<-as.numeric(dip_pred)
  diploid_predicted[,i]<-dip_pred
  
  
  dip_t<-paste(res[,hap1+2],res[,hap2+2],sep="")
  dip_t<- gsub("32","23",gsub("31","13",gsub("21","12",dip_t)))
  dip_t[dip_t == "11"] <- 1
  dip_t[dip_t == "12"] <- 2
  dip_t[dip_t == "13"] <- 3
  dip_t[dip_t == "22"] <- 4
  dip_t[dip_t == "23"] <- 5
  dip_t[dip_t == "33"] <- 6
  dip_t<-as.numeric(dip_t)
  diploid_true[,i]<-dip_t
}
long<-cbind.data.frame(snps,diploid_predicted)
res<-cbind.data.frame(snps,diploid_true)
#long<-long + 1
#colnames(long)<-c("chm","pos",ids)
corrs<-c(rep(0,20))
str(long)
str(res)
for (i in c(1:20)){
  j<-i+2
  corrs[i]<-cor.test(as.numeric(unlist(long[,..j])),as.numeric(unlist(res[,..j])),method="pearson")[[4]]
}
str(corrs)
corrs[is.na(corrs)]<-1.001
corr_R2<-corrs*corrs

print("writing corrs")
fwrite(as.list(corrs),args$out %&% "_dip_accuracy.txt",sep="\t")
print("writing corr R2")
str(corr_R2)
corr_R2<-as.list(corr_R2)
str(corr_R2)
fwrite(corr_R2,args$out %&% "_dip_accuracy.txt",sep="\t",append=T)