################
####Analyze family data by lagging versus non-lagging embryos
################

##################
###Outersect
##################

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

##################
###Family 1
##################

data<-read.csv(file="genotypes_CALL-III-20_mother_embryo_combined_individuals_mito.txt_focal_family1_transposed_reformat",sep="\t",head=TRUE)

data<-subset(data,data$chr=="ScyDAA6-1934-HRSCAF-2318")


focal<-c("E12_L002_R1_001.fastq","E4_L002_R1_001.fastq","E5_L002_R1_001.fastq","E6_L002_R1_001.fastq","E9_L002_R1_001.fastq","E.31.28_L004_R1_001.fastq","E.31.35_L004_R1_001.fastq","E.31.33_L004_R1_001.fastq")

focal_anc1<-data[focal]

non_focal<-outersect(colnames(data),c(focal,"chr","pos"))

comp_anc1<-data[non_focal]

mean_focal1<-rowMeans(focal_anc1,na.rm=TRUE)

mean_comp1<-rowMeans(comp_anc1,na.rm=TRUE)

##################
###Family 2
##################

data<-read.csv(file="genotypes_CALL-III-20_mother_embryo_combined_individuals_mito.txt_focal_family2",sep="\t",head=TRUE)

data<-subset(data,data$chr=="ScyDAA6-1934-HRSCAF-2318")

focal<-c("E12_L002_R1_001.fastq","E16_L002_R1_001.fastq","E11_L002_R1_001.fastq","E14_L002_R1_001.fastq","E21_L002_R1_001.fastq","E3_L002_R1_001.fastq","E8_L002_R1_001.fastq","E10_L002_R1_001.fastq","E18_L002_R1_001.fastq","E19_L002_R1_001.fastq","E.32.12_L004_R1_001.fastq")

focal_anc2<-data[focal]

non_focal<-outersect(colnames(data),c(focal,"chr","pos"))

comp_anc2<-data[non_focal]

mean_focal2<-rowMeans(focal_anc2,na.rm=TRUE)

mean_comp2<-rowMeans(comp_anc2,na.rm=TRUE)

##################
###Family 3
##################

data<-read.csv(file="genotypes_CALL-III-20_mother_embryo_combined_individuals_mito.txt_focal_family3",sep="\t",head=TRUE)

data<-subset(data,data$chr=="ScyDAA6-1934-HRSCAF-2318")

focal<-c("E11_L002_R1_001.fastq","E13_L002_R1_001.fastq","E14_L002_R1_001.fastq","E19_L002_R1_001.fastq","E6_L002_R1_001.fastq","E9_L002_R1_001.fastq","E.41.22_L004_R1_001.fastq")

focal_anc3<-data[focal]

non_focal<-outersect(colnames(data),c(focal,"chr","pos"))

comp_anc3<-data[non_focal]

mean_focal3<-rowMeans(focal_anc3,na.rm=TRUE)

mean_comp3<-rowMeans(comp_anc3,na.rm=TRUE)


##################
###Family 4
##################

data<-read.csv(file="genotypes_CALL-III-20_mother_embryo_combined_individuals_mito.txt_family4",sep="\t",head=TRUE)

data<-subset(data,data$chr=="ScyDAA6-1934-HRSCAF-2318")

focal<-c("E10_L002_R1_001.fastq","E11_L002_R1_001.fastq","E12_L002_R1_001.fastq","E21_L002_R1_001.fastq","E4_L002_R1_001.fastq","E20_L002_R1_001.fastq","E7_L002_R1_001.fastq","E.6.41_L005_R1_001.fastq","E.6.20_L005_R1_001.fastq")

focal_anc4<-data[focal]

non_focal<-outersect(colnames(data),c(focal,"chr","pos"))

comp_anc4<-data[non_focal]

mean_focal4<-rowMeans(focal_anc4,na.rm=TRUE)

mean_comp4<-rowMeans(comp_anc4,na.rm=TRUE)



##################
###Combine data and plot
##################

combined_geno_focal<-cbind(focal_anc1,focal_anc2,focal_anc3,focal_anc4)
combined_geno_comp<-cbind(comp_anc1,comp_anc2,comp_anc3,comp_anc4)

bir_focal<-{}
bir_comp<-{}
for(x in 1:length(combined_geno_focal[,1])){

focal<-length(subset(t(combined_geno_focal[x,]),t(combined_geno_focal[x,])==0))/length(t(combined_geno_focal[x,]))
comp<-length(subset(t(combined_geno_comp[x,]),t(combined_geno_comp[x,])==0))/length(t(combined_geno_comp[x,]))

bir_focal<-c(bir_focal,focal)
bir_comp<-c(bir_comp,comp)

}


plot(data$pos/1e6,bir_focal,type="l",xlim=c(1.9,2.2),xlab="Position",ylab="Proportion homozygous birchmanni",col="red",lwd=3,ylim=c(0,0.75))
points(data$pos/1e6,bir_comp,type="l",xlim=c(1.9,2.2),xlab="Position",ylab="Average ancestry",col="blue",lwd=3)
abline(v=2.066,lty=2)