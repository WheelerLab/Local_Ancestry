### Local ancestry acuracy ELAI
library(data.table)
library(argparse)
library(dplyr)
library(MOSAIC)
"%&%" = function(a,b) paste(a,b,sep="")

parser <- ArgumentParser()
parser$add_argument("--ps.21", help="resultant output of ELAI")
parser$add_argument("--haps.sample", help="admixed samle list as in MOSAIC intermediate files")
parser$add_argument("--nancestries", help="number of ancestries")
parser$add_argument("--result", help="results file output by adsimr")
parser$add_argument("--pos", help="snplist")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

#Read in snp ids
print("processing snp ids")
bim<-fread(args$pos, header = F, drop = "V1")
#dim(bim)
colnames(bim)<-c("pos","chm")
bim<-select(bim,chm,pos)
elai_out<-fread(args$ps.21,header = F)

#fread(args$ps.21, col.names = F)

#read in true ancestry and keep only snps found in bim file
print("reading in  true ancestry")
res<-as.data.frame(fread(args$result,header = T, showProgress = T))
res<-inner_join(bim,res,by=c("chm","pos"))

#create two empty data frames to fill with diploid ancestries
diploid_true<-data.frame(matrix(ncol = 20,nrow=50000))
local_df<-data.frame(matrix(ncol = 20,nrow=50000))
#convert haploid ancestries to diploid ancestries
if (args$nancestries == 3){
  ncols<-ncol(elai_out)
  random<-c(1,2,3)
  anc1_index<-seq(1, ncols, 3)
  anc2_index<-seq(2, ncols, 3)
  anc3_index<-seq(3, ncols, 3)
  for ( i in c(1:20)){ 
    cat( i, "/20\n")
    hap2<-i*2
    hap1<-hap2-1
    

    anc1_dosage<-round(elai_out[..i,..anc1_index],0)
    anc2_dosage<-round(elai_out[..i,..anc2_index],0)
    anc3_dosage<-round(elai_out[..i,..anc3_index],0)
    anc1_dosage[anc1_dosage==2]<-11
    anc2_dosage[anc2_dosage==2]<-22
    anc3_dosage[anc3_dosage==2]<-33
    anc1_dosage[anc1_dosage==1]<-1
    anc2_dosage[anc2_dosage==1]<-2
    anc3_dosage[anc3_dosage==1]<-3
    inidv<-cbind.data.frame(anc1_dosage,cbind.data.frame(anc2_dosage,anc3_dosage))
    indiv<-apply(indiv,1,sort,decreasing=F) %>% t()
    indiv_diploid<-indiv[,1] %&% indiv[,2] %&% indiv[,3]
    indiv_diploid<-as.numeric(indiv_diploid)
    indiv_diploid[indiv_diploid==0]<-sample(random,1)
    indiv_diploid[indiv_diploid == 11] <- 1
    indiv_diploid[indiv_diploid == 12] <- 2
    indiv_diploid[indiv_diploid == 13] <- 3
    indiv_diploid[indiv_diploid == 22] <- 4
    indiv_diploid[indiv_diploid == 23] <- 5
    indiv_diploid[indiv_diploid == 33] <- 6
    local_df[,i]<-indiv_diploid
    
    #repeat process for results df
    dip_t<-paste(res[,hap1+2],res[,hap2+2],sep="")
    dip_t<-gsub("32","23",gsub("31","13",gsub("21","12",dip_t)))
    dip_t<-gsub("32","23",gsub("31","13",gsub("21","12",dip_t)))
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
} else if (args$nancestries ==2){
  random<-c(1,2)
  ncols<-ncol(elai_out)
  anc1_index<-seq(1, ncols, 2)
  anc2_index<-seq(2, ncols, 2)
  anc1_population<-round(elai_out[,..anc1_index],0) %>% t()
  anc2_population<-round(elai_out[,..anc2_index],0) %>% t()
  str(anc1_population)
  str(anc2_population)
  for ( i in c(1:20)){ 
    cat( i, "/ 20\n")
    hap2<-i*2
    hap1<-hap2-1
    
    print("getting ancestral dosages")
    anc1_dosage<-anc1_population[,i]
    anc2_dosage<-anc2_population[,i]
    print("encoding")
    anc1_dosage[anc1_dosage==2]<-11
    anc2_dosage[anc2_dosage==2]<-22
    anc1_dosage[anc1_dosage==1]<-1
    anc2_dosage[anc2_dosage==1]<-2
    print("creating diploid encoding")
    indiv<-cbind.data.frame(anc1_dosage,anc2_dosage)
    indiv<-apply(indiv,1,sort,decreasing=F) %>% t()
    indiv_diploid<-indiv[,1] %&% indiv[,2]
    indiv_diploid<-as.numeric(indiv_diploid) 
    indiv_diploid[indiv_diploid == 11] <- 1
    indiv_diploid[indiv_diploid == 12] <- 2
    indiv_diploid[indiv_diploid == 22] <- 3
    local_df[,i]<-indiv_diploid
    print("creating dipoid true ancestry")
    dip_t<-paste(res[,hap1+2],res[,hap2+2],sep="")
    dip_t<-gsub("21","12",dip_t)
    dip_t[dip_t == "11"] <- 1
    dip_t[dip_t == "12"] <- 2
    dip_t[dip_t == "22"] <- 3
    diploid_true[,i]<-as.numeric(dip_t)
  }
}
## assign snp ids to MOSAIC
local_df<-as.data.frame(cbind.data.frame(bim,local_df))
#colnames(local_df)<-c("chm","pos",ids)
diploid_true<-as.data.frame(cbind.data.frame(bim,diploid_true))
##read in known ancestry, make sure it has only the right snp, and cols in the right order
#fwrite(local_df)

#sanity check
print("MOS dim:")
dim(local_df)
# str(local_df)
print("actual dim")
dim(diploid_true)
# str(diploid_true)
# fwrite(diploid_true,"~/software/Local_Ancestry/MOS_test_True.txt", sep = "\t")
# fwrite(local_df,"~/software/Local_Ancestry/MOS_test_estimated.txt", sep = "\t")
#begin corr tests per haplotype
corrs<-c(rep(0,20))
for (i in c(1:20)){
  corrs[i]<-cor.test(local_df[,(i+2)],diploid_true[,(i+2)],method="pearson")[[4]]
  # print(sd(local_df[,(i+2)]))
  # print(sd(diploid_true[,(i+2)]))
}
warnings()
corrs[is.na(corrs)]<-1.01
str(corrs)
#fwrite(local_df,"~/software/Local_Ancestry/Local_DF.txt",col.names = T,sep='\t' )
#process and write to file
corr_R2<-corrs*corrs

print("writing corrs")
fwrite(as.list(corrs),args$out %&% "_dip_accuracy.txt",sep="\t")
print("writing corr R2")
corr_R2<-as.list(corr_R2)
str(corr_R2)
fwrite(corr_R2,args$out %&% "_dip_accuracy.txt",sep="\t",append=T)