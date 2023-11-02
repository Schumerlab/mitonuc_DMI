library(tidyverse)
library(ggtext)

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

mitocopydata<-read.csv(file="mitoqpcr_ND1_Nup43_diluted.liver.DNA.csv") %>%
  mutate(copy_ratio = 2^delta_Ct)
mitocopydata$group <- factor(mitocopydata$group, levels = c("xb","f1","xm"))
summary(aov(mitocopydata$delta_Ct~as.factor(mitocopydata$group)))
TukeyHSD(aov(mitocopydata$delta_Ct~as.factor(mitocopydata$group)))

mitocopies <- ggplot(mitocopydata,aes(x=group, y=delta_Ct)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  #scale_y_reverse() +
  #scale_y_continuous(trans = "log2") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,hetcol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,hetcol,malcol), width=0.15,size=0.8)  + 
  #labs(x = "Genotype", y = expression("qPCR mtDNA : nDNA ("~C[t]^n~-C[t]^mt~")")) +
  labs(x = "Genotype", y = expression("qPCR mtDNA : nDNA ("~-Delta*C[t]~")")) +
  scale_x_discrete(labels = c("*X. birchmanni*", "hybrid", "*X. malinche*"))
mitocopies
