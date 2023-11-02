a=0.5
mm=a^2
mb=2*a*(1-a)
bb=(1-a)^2

simulation_results<-{}

for(x in 1:500000){

s=runif(1,0,1)
h=runif(1,0,1)

N=10000

fmm=mm*N
fmb=mb*N*(1-h*s)
fbb=bb*N*(1-s)

N_total=fmm+fmb+fbb

ratio_bb=fbb/N_total
ratio_mb=fmb/N_total
ratio_mm=fmm/N_total

pop<-c(rep(2,fmm),rep(1,fmb),rep(0,fbb))

N_indiv=937

sim_pop=sample(pop,N_indiv)

obb<-length(subset(sim_pop,sim_pop==0))
obm<-length(subset(sim_pop,sim_pop==1))

a<- (obb + 0.5*obm)/N_indiv

simulation_results<-rbind(simulation_results,cbind(s,h,obb,obm,a))

}

accepted<-subset(simulation_results,simulation_results[,5] <(0.31+0.015) & simulation_results[,5] > (0.31-0.015) & simulation_results[,3] <=1 & simulation_results[,4] <=611 & simulation_results[,4]>=593)
#all +/- binomial standard error

write.table(accepted,file="~/Box/Schumer_lab_resources/Project_files/Chromosome13_incompatibility/Data/ABC_simulations_accepted_F2s_March2021.txt",sep="\t",row.names=FALSE,quote=FALSE)