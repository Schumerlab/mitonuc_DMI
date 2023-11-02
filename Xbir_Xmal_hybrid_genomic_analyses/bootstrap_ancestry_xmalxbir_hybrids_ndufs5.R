tlmc<-read.csv(file="TLMC_50kb_window_ancestry_2017",sep="\t",head=FALSE)
acua<-read.csv(file="ACUA_50kb_window_ancestry_2018",sep="\t",head=FALSE)
agcz<-read.csv(file="AGCZ_50kb_window_ancestry_multiyear",sep="\t",head=FALSE)

###reset tlmc so minor parent ancestry is polarized the same way in each population
tlmc$V4<- 1-tlmc$V4

###set observed minor parent ancestry at ndufs5
obs1=0.050
obs2=0.022
obs3=0.021


### First test by simple random selection of ancestry windows
tlmc<-na.omit(tlmc)
acua<-na.omit(acua)
agcz<-na.omit(acua)

simulations<-{}

for(x in 1:10000){

rand_tlmc<-sample(tlmc$V4,1)
rand_acua<-sample(acua$V4,1)
rand_agcz<-sample(agcz$V4,1)

simulations<-rbind(simulations,cbind(rand_tlmc,rand_acua,rand_agcz))

}

subset(simulations,simulations[,1]<=obs1  & simulations[,2]<=obs2 & simulations[,3]<=obs3)



#### Now with a sliding window approach (maintaining linkage structure)

tlmc<-read.csv(file="TLMC_50kb_window_ancestry_2017",sep="\t",head=FALSE)
acua<-read.csv(file="ACUA_50kb_window_ancestry_2018",sep="\t",head=FALSE)
agcz<-read.csv(file="AGCZ_50kb_window_ancestry_multiyear",sep="\t",head=FALSE)

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

megaframe <- agcz %>%
  rename(`chr` = V1, `pos` = V2, `end` = V3, `anc` = V4, `whatever` = V5, `blah` = V6) %>%
  mutate(tlmc_chr = tlmc$V1, tlmc_pos = tlmc$V2, tlmc_anc = tlmc$V4, acua_chr = acua$V1, acua_pos = acua$V2, acua_anc = acua$V4)

megaframe_null <- megaframe

filter(megaframe, anc < 0.021, acua_anc < 0.022, tlmc_anc > 0.94)

step1 = 100
step2 = 200

nulls<-{}
comparisons <- 0

total <- nrow(megaframe)

megaframe_null <- megaframe_null %>%
  mutate(tlmc_chr = shifter(tlmc_chr, n = step1), tlmc_pos = shifter(tlmc_pos, n = step1), tlmc_anc = shifter(tlmc_anc, n = step1),
         acua_chr = shifter(acua_chr, n = step2), acua_pos = shifter(acua_pos, n = step2), acua_anc = shifter(acua_anc, n = step2))

nulls <- rbind(nulls, filter(megaframe_null, anc < 0.021, acua_anc < 0.022, tlmc_anc > 0.94))
comparisons <- comparisons + total

while (any(megaframe_null[1,] != megaframe[1,])) {
  megaframe_null <- megaframe_null %>%
    mutate(tlmc_chr = shifter(tlmc_chr, n = step1), tlmc_pos = shifter(tlmc_pos, n = step1), tlmc_anc = shifter(tlmc_anc, n = step1),
           acua_chr = shifter(acua_chr, n = step2), acua_pos = shifter(acua_pos, n = step2), acua_anc = shifter(acua_anc, n = step2))
  
  nulls <- rbind(nulls, filter(megaframe_null, anc < 0.021, acua_anc < 0.022, tlmc_anc > 0.94))
  comparisons <- comparisons + total
}

nrow(nulls)/comparisons