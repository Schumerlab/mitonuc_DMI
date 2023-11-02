#several simple power simulations for nuclear-nuclear incompatibilities

#####symmetric
index<-read.csv(file="./Data/hybrid_index_allindivs_CALL_Nov2021",sep="\t",head=TRUE)
index<-na.omit(index)

s=0.998
num_sims=2000
num_indiv=339

power_results<-{}

for(k in 1:1000){
index_total<-rep(index$hybrid_index,ceiling(num_sims/length(index$hybrid_index)))
nuclear_sim1<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)
nuclear_sim2<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)


geno_post_selection<-{}

for(x in 1:length(nuclear_sim1)){

nuclear_sim1_focal=nuclear_sim1[x]
nuclear_sim2_focal=nuclear_sim2[x]

if(((nuclear_sim1_focal==0) & (nuclear_sim2_focal==2)) | ((nuclear_sim1_focal==2) & (nuclear_sim2_focal==0))){

survive=rbinom(1,1,1-s)
if(survive==1){geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))}
} else{
geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))
}

}


geno_sample<-geno_post_selection[sample(nrow(geno_post_selection),num_indiv),]

power_results<-c(power_results,pcor(geno_sample)$p.value[,1][2])

}

length(subset(power_results,power_results<4e-5))

#####symmetric dominant
index<-read.csv(file="./Data/hybrid_index_allindivs_CALL_Nov2021",sep="\t",head=TRUE)
index<-na.omit(index)

s=0.998
num_sims=2000
num_indiv=339

power_results<-{}

for(k in 1:1000){
index_total<-rep(index$hybrid_index,ceiling(num_sims/length(index$hybrid_index)))
nuclear_sim1<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)
nuclear_sim2<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)


geno_post_selection<-{}

for(x in 1:length(nuclear_sim1)){

nuclear_sim1_focal=nuclear_sim1[x]
nuclear_sim2_focal=nuclear_sim2[x]

if(((nuclear_sim1_focal==0) & (nuclear_sim2_focal==2)) | ((nuclear_sim1_focal==2) & (nuclear_sim2_focal==0)) | ((nuclear_sim1_focal==1) & (nuclear_sim2_focal==0)) | ((nuclear_sim1_focal==2) & (nuclear_sim2_focal==1)) | ((nuclear_sim1_focal==1) & (nuclear_sim2_focal==1))){

survive=rbinom(1,1,1-s)
if(survive==1){geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))}
} else{
geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))
}

}


geno_sample<-geno_post_selection[sample(nrow(geno_post_selection),num_indiv),]

power_results<-c(power_results,pcor(geno_sample)$p.value[,1][2])

}

length(subset(power_results,power_results<4e-5))

#####asymmetric dominant
s=0.998
num_sims=2000
num_indiv=339

power_results<-{}

for(k in 1:1000){
index_total<-rep(index$hybrid_index,ceiling(num_sims/length(index$hybrid_index)))
nuclear_sim1<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)
nuclear_sim2<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)


geno_post_selection<-{}

for(x in 1:length(nuclear_sim1)){

nuclear_sim1_focal=nuclear_sim1[x]
nuclear_sim2_focal=nuclear_sim2[x]

if(((nuclear_sim1_focal==2) & (nuclear_sim2_focal==0)) | ((nuclear_sim1_focal==1) & (nuclear_sim2_focal==0)) | ((nuclear_sim1_focal==1) & (nuclear_sim2_focal==1))){

survive=rbinom(1,1,1-s)
if(survive==1){geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))}
} else{
geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))
}

}


geno_sample<-geno_post_selection[sample(nrow(geno_post_selection),num_indiv),]

power_results<-c(power_results,pcor(geno_sample)$p.value[,1][2])

}

length(subset(power_results,power_results<4e-5))

#####asymmetric recessive
s=0.998
num_sims=2000
num_indiv=339

power_results<-{}

for(k in 1:1000){
index_total<-rep(index$hybrid_index,ceiling(num_sims/length(index$hybrid_index)))
nuclear_sim1<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)
nuclear_sim2<-rbinom(num_sims,1,index$hybrid_index)+rbinom(num_sims,1,index$hybrid_index)


geno_post_selection<-{}

for(x in 1:length(nuclear_sim1)){

nuclear_sim1_focal=nuclear_sim1[x]
nuclear_sim2_focal=nuclear_sim2[x]

if(((nuclear_sim1_focal==2) & (nuclear_sim2_focal==0))){

survive=rbinom(1,1,1-s)
if(survive==1){geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))}
} else{
geno_post_selection<-rbind(geno_post_selection,cbind(nuclear_sim1_focal,nuclear_sim2_focal,index_total[x]))
}

}


geno_sample<-geno_post_selection[sample(nrow(geno_post_selection),num_indiv),]

power_results<-c(power_results,pcor(geno_sample)$p.value[,1][2])

}

length(subset(power_results,power_results<4e-5))