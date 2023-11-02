#several simple power simulations for revisions:

library("ppcor")

index<-read.csv(file="hybrid_index_allindivs_CALL_Nov2021",sep="\t",head=TRUE)

#################
#####asymmetric
#################
s=0.9
mito_sel=1
nuclear_sel=0
num_sims=2000
num_indiv=330

power_results<-{}

for(k in 1:1000){
index_total<-rep(index$hybrid_index,ceiling(num_sims/length(index$hybrid_index)))
mito_sim<-rbinom(num_sims,1,index$hybrid_index)
nuclear_sim<-rbinom(num_sims,1,index$hybrid_index)+rbinom(10000,1,index$hybrid_index)

geno_post_selection<-{}

for(x in 1:length(mito_sim)){

mito_focal=mito_sim[x]
nuclear_focal=nuclear_sim[x]

if((mito_focal==mito_sel) & (nuclear_focal==nuclear_sel)){

survive=rbinom(1,1,1-s)
if(survive==1){geno_post_selection<-rbind(geno_post_selection,cbind(mito_focal,nuclear_focal,index_total[x]))}
} else{
geno_post_selection<-rbind(geno_post_selection,cbind(mito_focal,nuclear_focal,index_total[x]))
}

}


geno_sample<-geno_post_selection[sample(nrow(geno_post_selection),num_indiv),]

power_results<-c(power_results,pcor(geno_sample)$p.value[,1][2])

}

####################
####symmetric
####################

index<-read.csv(file="./Chromosome13_incompatibility/Data/hybrid_index_allindivs_CALL_Nov2021",sep="\t",head=TRUE)

#####asymmetric
s1=0.998
s2=0.17
num_sims=2000
num_indiv=330


power_results_sym<-{}

for(k in 1:1000){
index_total<-rep(index$hybrid_index,ceiling(num_sims/length(index$hybrid_index)))
mito_sim<-rbinom(num_sims,1,index$hybrid_index)
nuclear_sim<-rbinom(num_sims,1,index$hybrid_index)+rbinom(10000,1,index$hybrid_index)

geno_post_selection<-{}

for(x in 1:length(mito_sim)){

mito_focal=mito_sim[x]
nuclear_focal=nuclear_sim[x]

if(((mito_focal==0) & (nuclear_focal==2)) | ((mito_focal==1) & (nuclear_focal==0))){

if((mito_focal==1) & (nuclear_focal==0)){survive=rbinom(1,1,1-s1)}
if((mito_focal==0) & (nuclear_focal==2)){survive=rbinom(1,1,1-s2)}
if(survive==1){geno_post_selection<-rbind(geno_post_selection,cbind(mito_focal,nuclear_focal,index_total[x]))}
} else{
geno_post_selection<-rbind(geno_post_selection,cbind(mito_focal,nuclear_focal,index_total[x]))
} 

}

geno_sample<-geno_post_selection[sample(nrow(geno_post_selection),num_indiv),]

power_results_sym<-c(power_results_sym,pcor(geno_sample)$p.value[,1][2])

}
