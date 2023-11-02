library(tidyverse)
library(cowplot)

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

slides <- read.csv("input_files/F2_NDUFA13_histology_slides_analysis_cleaned.csv") %>%
  mutate(incompatible = ndufa13 == 0,
         compatible = ndufa13 == 2,
         total_area_um = total_area * (.5049) ^ 2,
         total_area_mm = total_area_um / 1000000)

maxes <- filter(slides, max_area_bool == T)

t.test(filter(maxes, incompatible == F, type == "atrium", damaged == F)$total_area, filter(maxes, incompatible == T, type == "atrium", damaged == F)$total_area)  
t.test(filter(maxes, incompatible == F, type == "atrium", damaged == F)$mc_percent_quad, filter(maxes, incompatible == T, type == "atrium", damaged == F)$mc_percent_quad)  

atr_totalarea <- ggplot(filter(maxes, type == "atrium", damaged == F),aes(x=compatible, y=total_area_mm)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  #ggtitle("Atrium") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,malcol), width=0.15,size=0.8)  + 
  labs(x = "", y = expression("Cross-Sectional Area ("*mm^2*")")) +
  scale_x_discrete(labels = c("", ""))
atr_totalarea

vnt_totalarea <- ggplot(filter(maxes, type == "ventricle", damaged == F),aes(x=compatible, y=total_area_mm)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  #ggtitle("Ventricle") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,malcol), width=0.15,size=0.8)  + 
  labs(x = expression(italic("ndufa13")), y = expression("Cross-Sectional Area ("*mm^2*")")) +
  scale_x_discrete(labels = c("*X. birchmanni*", "*X. malinche*"))
vnt_totalarea

atr_muscle <- ggplot(filter(maxes, type == "atrium", damaged == F),aes(x=compatible, y=mc_percent_quad)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,malcol), width=0.15,size=0.8)  + 
  labs(x = "", y = expression("Area of Cardiomyocytes (%)")) +
  scale_x_discrete(labels = c("", ""))
atr_muscle

vnt_muscle <- ggplot(filter(maxes, type == "ventricle", damaged == F),aes(x=incompatible, y=mc_percent_quad)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,malcol), width=0.15,size=0.8)  + 
  #labs(x = "Genotype", y = expression("qPCR mtDNA : nDNA ("~C[t]^n~-C[t]^mt~")")) +
  labs(x = expression(italic("ndufa13")), y = expression("Area of Cardiomyocytes (%)")) +
  scale_x_discrete(labels = c("*X. birchmanni*", "*X. malinche*"))
vnt_muscle


atr_blood <- ggplot(filter(maxes, type == "atrium", damaged == F),aes(x=compatible, y=bc_percent_quad)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,malcol), width=0.15,size=0.8)  + 
  labs(x = expression(""), y = expression("Area of Blood Cells (%)")) +
  scale_x_discrete(labels = c("", ""))
atr_blood

vnt_blood <- ggplot(filter(maxes, type == "ventricle", damaged == F),aes(x=incompatible, y=bc_percent_quad)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.2,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,malcol), width=0.15,size=0.8)  + 
  labs(x = expression(italic("ndufa13")), y = expression("Area of Blood Cells (%)")) +
  scale_x_discrete(labels = c("*X. birchmanni*", "*X. malinche*"))
vnt_blood

pic_standin <- ggplot() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_blank()) +
  geom_blank()

topplot <- plot_grid(atr_totalarea, atr_muscle, atr_blood, vnt_totalarea, vnt_muscle, vnt_blood,
                       labels = c("a", "b", "c", "d", "e", "f"), ncol = 3, align = "v")
bottomplot <- plot_grid(pic_standin, pic_standin, labels = c("G", "H"), nrow = 1)
wholeplot <- plot_grid(topplot, bottomplot, nrow = 2, rel_heights = c(2,1))
wholeplot

ggsave("output_files/histology_cowplot.pdf",
       wholeplot, width = 9, height = 9, units = "in")


mitocopies