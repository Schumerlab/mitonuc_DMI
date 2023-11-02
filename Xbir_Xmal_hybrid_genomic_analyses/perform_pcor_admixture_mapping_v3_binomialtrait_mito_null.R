######perform admixture mapping given a genotypes file, a file with hybrid index (expects hybrid index to be in the second column, a phenotypes file, and a column number for the phenotype of interest
#####NOTE: individual order *must* be identical across files!!!

args <- commandArgs(TRUE)
if(length(args)<5){
stop("usage is Rscript perform_pcor_admixture_mapping_v3_binomialtrait_mito_null.R genotypes_file hybrid_index_file phenotypes_file focal_column_number name_tag")
}

library("ppcor")

infile <- as.character(args[1])

data<-read.csv(file=infile,head=T,as.is=T,sep="\t")

hybrid_index<-as.character(args[2])

pheno<-as.character(args[3])

pheno_column<-as.numeric(args[4])

phenotypes<-read.csv(file=pheno,sep="\t",head=TRUE)

index<-read.csv(file=hybrid_index,sep="\t",head=TRUE)

tag<-as.character(args[5])

for(k in 1:500){
out<-paste(infile,"_results_ppcor_v2_null_",tag,"_",k,sep="")
file.remove(out)

names<-colnames(data)

track=0
for (x in 2:length(data[1,])){

dat<-cbind(data[,x],phenotypes[,pheno_column],index$hybrid_index)
dat<-na.omit(dat)

rand_pheno<-rbinom(length(dat[,3]),1,dat[,3])
dat[,2]<-rand_pheno

if(length(unique(dat[,2]))>1){

pval<-pcor(dat)$p.value[2]
estimate<-pcor(dat)$estimate[2]

results<-cbind(names[x],pval,estimate,length(dat[,1]))
if(track==0){
write.table(results,file=out,append=TRUE,col.names=c("chrom.marker","pval","estimate","num_ind"),row.names=F,sep="\t",quote=FALSE)
track=1
} else{
write.table(cbind(names[x],pval,estimate,length(dat[,1])),file=out,append=TRUE,col.names=FALSE,row.names=F,sep="\t",quote=FALSE)

}

}#if half the data

}#for all lines

warnings()

}