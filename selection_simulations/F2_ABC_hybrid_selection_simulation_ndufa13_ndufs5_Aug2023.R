a=0.5
mm=a^2
mb=2*a*(1-a)
bb=(1-a)^2

simulation_results<-{}

posterior1<-read.csv(file="ABC_simulations_accepted_F2s_March2021.txt",sep="\t",head=TRUE)
posterior2<-read.csv(file="ABC_simulations_accepted_F2s_March2021_ndufa13.txt",sep="\t",head=TRUE)

for(x in 1:10000){

N=10000

#simulate selection on locus1
r1<-posterior1[sample(nrow(posterior1),1),]
s1<-r1[,1]
h1<-r1[,2]
fmm1=mm*N
fmb1=mb*N*(1-h1*s1)
fbb1=bb*N*(1-s1)

N_total1=fmm1+fmb1+fbb1

ratio_bb1=fbb1/N_total1
ratio_mb1=fmb1/N_total1
ratio_mm1=fmm1/N_total1

#populate locus1
pop1<-c(rep(2,fmm1),rep(1,fmb1),rep(0,fbb1))

#simulate selection on locus2
r2<-posterior2[sample(nrow(posterior2),1),]
s2<-r2[,1]
h2<-r2[,2]
fmm2=mm*N
fmb2=mb*N*(1-h2*s2)
fbb2=bb*N*(1-s2)

N_total2=fmm2+fmb2+fbb2

ratio_bb2=fbb2/N_total2
ratio_mb2=fmb2/N_total2
ratio_mm2=fmm2/N_total2

#populate locus2
pop2<-c(rep(2,fmm2),rep(1,fmb2),rep(0,fbb2))

N_indiv=1010 #num with data at both loci

sim_pop1=sample(pop1,N_indiv)
sim_pop2=sample(pop2,N_indiv)

obb1<-length(subset(sim_pop1,sim_pop1==0))
obm1<-length(subset(sim_pop1,sim_pop1==1))

obb2<-length(subset(sim_pop2,sim_pop2==0))
obm2<-length(subset(sim_pop2,sim_pop2==1))

a1<- (obb1 + 0.5*obm1)/N_indiv
a2<- (obb2 + 0.5*obm2)/N_indiv

two_locus1<-length(subset(sim_pop2,sim_pop2==0 & sim_pop1==2))
two_locus2<-length(subset(sim_pop2,sim_pop2==0 & sim_pop1!=2))

simulation_results<-rbind(simulation_results,cbind(s1,h1,s2,h2,obb1,obm1,obb2,obm2,a1,a2,two_locus1,two_locus2))

}
