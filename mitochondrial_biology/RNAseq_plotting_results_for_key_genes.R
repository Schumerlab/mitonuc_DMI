##############
###RNAseq results for key genes, data from Payne et al. Molecular Ecology 2022
##############

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
 if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
 stop("vectors must be same length")
 	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
 	}

library("plotrix")

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

#########
###ndufs5
#########

data<-read.csv(file="RNAseq_data_Payne_etal_mitoDMI.csv",head=TRUE)

f1_ndufs5<-subset(data$Normalized_counts,data$gene=="ndufs5" & data$genotype=="f1")
xbir_ndufs5<-subset(data$Normalized_counts,data$gene=="ndufs5" & data$genotype=="xbir")
xmal_ndufs5<-subset(data$Normalized_counts,data$gene=="ndufs5" & data$genotype=="xmal")

plot(1:3,c(mean(xmal_ndufs5),mean(f1_ndufs5),mean(xbir_ndufs5)),col=c(malcol,hetcol,bircol),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=20,cex=2,xlim=c(0.5,3.5),ylim=c(100,350))

error.bar(1:3,c(mean(xmal_ndufs5),mean(f1_ndufs5),mean(xbir_ndufs5)),2*c(std.error(xmal_ndufs5),std.error(f1_ndufs5),std.error(xbir_ndufs5)),col=c(malcol,hetcol,bircol),lwd=2)

malcol=rgb(0/255,0/255,139/255,alpha=0.3)
hetcol=rgb(65/255,105/255,225/255,0.6)
bircol=rgb(255/255,0/255,0/255,alpha=0.6)

BB<-xbir_ndufs5
MB<-f1_ndufs5
MM<-xmal_ndufs5

noise<-runif(length(BB),0.25,0.35)
points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)

noise<-runif(length(MB),0.25,0.35)
points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)

noise<-runif(length(MM),0.25,0.35)
points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)

mtext(c("MM","MB","BB"),at=1:3,side=1)

##########
###ndufa13
###########

f1_ndufa13<-subset(data$Normalized_counts,data$gene=="ndufa13" & data$genotype=="f1")
xbir_ndufa13<-subset(data$Normalized_counts,data$gene=="ndufa13" & data$genotype=="xbir")
xmal_ndufa13<-subset(data$Normalized_counts,data$gene=="ndufa13" & data$genotype=="xmal")

plot(1:3,c(mean(xmal_ndufa13),mean(f1_ndufa13),mean(xbir_ndufa13)),col=c(malcol,hetcol,bircol),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=20,cex=2,xlim=c(0.5,3.5),ylim=c(100,450))

error.bar(1:3,c(mean(xmal_ndufa13),mean(f1_ndufa13),mean(xbir_ndufa13)),2*c(std.error(xmal_ndufa13),std.error(f1_ndufa13),std.error(xbir_ndufa13)),col=c(malcol,hetcol,bircol),lwd=2)

malcol=rgb(0/255,0/255,139/255,alpha=0.3)
hetcol=rgb(65/255,105/255,225/255,0.6)
bircol=rgb(255/255,0/255,0/255,alpha=0.6)

BB<-xbir_ndufa13
MB<-f1_ndufa13
MM<-xmal_ndufa13

noise<-runif(length(BB),0.25,0.35)
points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)

noise<-runif(length(MB),0.25,0.35)
points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)

noise<-runif(length(MM),0.25,0.35)
points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)

mtext(c("MM","MB","BB"),at=1:3,side=1)

#########
##nd2
#########

f1_nd2<-subset(data$Normalized_counts,data$gene=="nd2" & data$genotype=="f1")
xbir_nd2<-subset(data$Normalized_counts,data$gene=="nd2" & data$genotype=="xbir")
xmal_nd2<-subset(data$Normalized_counts,data$gene=="nd2" & data$genotype=="xmal")

plot(1:3,c(mean(xmal_nd2),mean(f1_nd2),mean(xbir_nd2)),col=c(malcol,hetcol,bircol),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=20,cex=2,xlim=c(0.5,3.5),ylim=c(-10,8000))

error.bar(1:3,c(mean(xmal_nd2),mean(f1_nd2),mean(xbir_nd2)),2*c(std.error(xmal_nd2),std.error(f1_nd2),std.error(xbir_nd2)),col=c(malcol,hetcol,bircol),lwd=2)

malcol=rgb(0/255,0/255,139/255,alpha=0.3)
hetcol=rgb(65/255,105/255,225/255,0.6)
bircol=rgb(255/255,0/255,0/255,alpha=0.6)

BB<-xbir_nd2
MB<-f1_nd2
MM<-xmal_nd2

noise<-runif(length(BB),0.25,0.35)
points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)

noise<-runif(length(MB),0.25,0.35)
points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)

noise<-runif(length(MM),0.25,0.35)
points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)

mtext(c("MM","MB","BB"),at=1:3,side=1)

########
##nd6
#########

f1_nd6<-subset(data$Normalized_counts,data$gene=="nd6" & data$genotype=="f1")
xbir_nd6<-subset(data$Normalized_counts,data$gene=="nd6" & data$genotype=="xbir")
xmal_nd6<-subset(data$Normalized_counts,data$gene=="nd6" & data$genotype=="xmal")


plot(1:3,c(mean(xmal_nd6),mean(f1_nd6),mean(xbir_nd6)),col=c(malcol,hetcol,bircol),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=20,cex=2,xlim=c(0.5,3.5),ylim=c(0,3500))

error.bar(1:3,c(mean(xmal_nd6),mean(f1_nd6),mean(xbir_nd6)),2*c(std.error(xmal_nd6),std.error(f1_nd6),std.error(xbir_nd6)),col=c(malcol,hetcol,bircol),lwd=2)

malcol=rgb(0/255,0/255,139/255,alpha=0.3)
hetcol=rgb(65/255,105/255,225/255,0.6)
bircol=rgb(255/255,0/255,0/255,alpha=0.6)

BB<-xbir_nd6
MB<-f1_nd6
MM<-xmal_nd6

noise<-runif(length(BB),0.25,0.35)
points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)

noise<-runif(length(MB),0.25,0.35)
points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)

noise<-runif(length(MM),0.25,0.35)
points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)

mtext(c("MM","MB","BB"),at=1:3,side=1)
