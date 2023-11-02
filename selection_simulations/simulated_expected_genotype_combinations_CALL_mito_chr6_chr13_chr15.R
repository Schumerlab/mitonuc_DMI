#####simulate expections chr13 and chr6

data<-read.csv(file="genotypes_allchrs_Nov2021_CALL_v2.txt_focal_chr6_chr13",sep="\t",head=TRUE)
index<-read.csv(file="hybrid_index_allindivs_CALL_v2_Nov2021",sep="\t",head=TRUE)

x<-cbind(data[,2],data[,3],data[,4],index$hybrid_index)
colnames(x)<-c("chr6","chr13","mito","index")
x<-na.omit(x)

positives1<-{}
positives2<-{}
for(k in 1:10000){
geno<- rbinom(length(x[,4]),1,x[,4]) + rbinom(length(x[,4]),1,x[,4])
mito<- x[,3]
positives1<-c(positives1,length(subset(geno,geno==0 & mito==2)))
positives2<-c(positives2,length(subset(geno,geno==2 & mito==0)))
}


####simulate expectations chr 15
data<-read.csv(file="genotypes_allchrs_Nov2021_CALL_v2.txt_focal_chr15_chr16",sep="\t",head=TRUE)
index<-read.csv(file="hybrid_index_allindivs_CALL_v2_Nov2021",sep="\t",head=TRUE)

x<-cbind(data[,2],data[,4],index$hybrid_index)
colnames(x)<-c("chr15","mito","index")
x<-na.omit(x)


positives1<-{}
positives2<-{}
for(k in 1:10000){
geno<- rbinom(length(x[,3]),1,x[,3]) + rbinom(length(x[,3]),1,x[,3])
mito<- x[,2]
positives1<-c(positives1,length(subset(geno,geno==0 & mito==2)))
positives2<-c(positives2,length(subset(geno,geno==2 & mito==0)))
}
