library(tidyverse)

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

#ndufs5 ASE
a1<-c(50,51,26)
a2<-c(56,61,18)
ase <- a1/(a1+a2)
data_rna <- data.frame(intervals = NA, peptide = "WLLPQSGEQPYK", percent_malinche = ase, indiv = sapply(1:length(ase), function(x) paste("ind",x, sep = "")), isotope = "RNA")
datasumm <- group_by(data_rna, isotope) %>%
  summarize(mean = mean(percent_malinche),
            stdev = sd(percent_malinche),
            se = stdev/sqrt(n()))

ase_plot <- ggplot(filter(datasumm, isotope %in% c("RNA")), aes(x = isotope, y = mean, fill = isotope, color = isotope)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x.bottom = element_text(color = "black", size = 11)) +
  theme(panel.grid = element_blank()) +
  geom_point(stat = "identity") +
  geom_jitter(data = filter(data_rna, isotope %in% c("RNA")), aes(y = percent_malinche), width = 0.4, color = "grey") +
  lims(y = c(0.35, 0.65)) +
  labs(x = "RNA", y = expression(paste("Proportion ",italic("X. malinche")))) +
  scale_color_manual(values = c(hetcol, rgb(50/255,50/255,50/255), rgb(150/255,150/255,10/255))) +
  scale_x_discrete(labels = expression(italic("ndufs5"))) +
  geom_hline(aes(yintercept = 0.5), lty = 2, color = "grey") + 
  geom_errorbar(aes(ymin = mean - 2 * se, ymax = mean + 2*se), width=0.15,size=0.8)
ase_plot


#ndufa13
a1<-c(126,96,48)
a2<-c(95,81,49)
a2/(a1+a2)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

plot(1:3,a1/(a1+a2),xlab="Individual",ylab="Proportion X. birchmanni allele",col="blue",ylim=c(0.25,0.75),pch=20,cex=2,xaxt="n",xlim=c(0.75,3.25))
ratio<-a1/(a1+a2)
se<-sqrt(ratio*(1-ratio)/(a1+a2))
abline(h=0.5,lty=2,col="lightgray",lwd=2)
error.bar(1:3,a1/(a1+a2),2*se,lwd=2,col="blue")
