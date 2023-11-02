stac<-read.csv(file="average_ancestry_by_site_STAC.txt_ancestry_cM_windows_0.1cM_xbirchmanni10x_WG",sep="\t",head=FALSE)

huex<-read.csv(file="average_ancestry_by_site_HUEX.txt_ancestry_cM_windows_0.1cM_xbirchmanni10x_WG",sep="\t",head=FALSE)

obs1=0.969
obs2=0.975

stac<-na.omit(stac)
huex<-na.omit(huex)

simulations<-{}

for(x in 1:10000){

rand_stac<-sample(stac$V4,1)
rand_huex<-sample(huex$V4,1)

simulations<-rbind(simulations,cbind(rand_stac,rand_huex))

}

subset(simulations,simulations[,1]>obs1  & simulations[,2]>obs2)