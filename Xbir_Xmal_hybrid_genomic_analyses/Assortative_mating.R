index<-read.csv(file="CALL_assortative_mating_analysis.csv",sep=",",head=TRUE)

pop_index<-read.csv(file="hybrid_index_CALL_allmitos",sep="\t",head=TRUE)

sim<-{}
for(x in 1:500){

maternal<-sample(index$mom_index,1)
mate<-sample(pop_index$hybrid_index,1)
off<-rnorm(1,mean=(maternal+mate)/2,sd=0.045) 

sim<-c(sim,maternal-off)

}