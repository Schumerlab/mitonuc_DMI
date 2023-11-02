library(tidyverse)

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

datum <- read.csv("input_files/FCF_data_051021_wdate.csv")
datum_cleaned <- datum %>%
  mutate(FCF = recode(FCF, CI = "CI Flux Control", 
                      CI_CII_eff = "CI + CII Efficiency", 
                      CI_eff = "CI Efficiency",
                      CII_A = "CII-A Flux Control",
                      CII_B = "CII-B Flux Control",
                      ETS_capacity = "ETS Capacity",
                      CIV = "Excess CIV Capacity",
                      Max_prot_norm_resp = "Maximum Respiration",
                      Leak_prot_norm_resp = "LEAK State Respiration"),
         Sex = recode(sex, male = "Male", female = "Female"),
         Genotype = recode_factor(species, birchmanni = "bir", `F1-hybrid-reverse` = "mt_bir", `F1-hybrid` = "mt_mal", malinche = "mal", .ordered = T)
  )
allplots <- ggplot(data=datum_cleaned,aes(x=Genotype, y=FCF_value)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_jitter(width=0.1,aes(color=Sex)) + 
  stat_summary(fun = mean, geom="point", size=1, color="black") + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color="black", width=0.1)  + 
  labs(x = "Genotype", y = "FCF Value (Mean ± 2 SE)") +
  facet_wrap(~FCF, scales="free")
allplots
ggsave("output_files/all_FCF_data_wleak.pdf",
       allplots, width = 8, height = 7, units = "in")

########
# Complex I Efficiency
########
cieff <- ggplot(data=filter(datum, FCF == "CI_eff", species != 'F1-hybrid-reverse'),aes(x=species, y=FCF_value)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.1,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,hetcol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,hetcol,malcol), width=0.15,size=0.8)  + 
  labs(x = "Genotype", y = "Complex I Efficiency") +
  scale_x_discrete(labels = c("*X. birchmanni*", "hybrid", "*X. malinche*"))
cieff
ggsave("output_files/complexI_efficiency_BMM_version.pdf",
       cieff, width = 3.75, height = 3, units = "in")  

######
# Time to peak post-ADP
######
data<-read.csv(file="input_files/O2K_time_to_peak.csv")
timetoadp_plot <- ggplot(data, aes(x = species, y = adp_time)) +
  theme_bw() +
  ylim(0,3000) +
  labs(x = "Genotype", y = "Seconds to Peak CI") +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.1,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,hetcol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,hetcol,malcol), width=0.15,size=0.8) +
  scale_x_discrete(labels = c("*X. birchmanni*", "hybrid", "*X. malinche*"))
timetoadp_plot
ggsave("output_files/time_to_adp_peak_v2.pdf",
       timetoadp_plot, width = 3.5, height = 3, units = "in")

########
# Maximum Respiration (Fig. 3E)
########

maxres<- ggplot(data=filter(datum, FCF == "Max_prot_norm_resp", species != 'F1-hybrid-reverse'),aes(x=species, y=FCF_value)) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.1,col="gray") +
  stat_summary(fun = mean, geom="point", size=5, shape=20, color=c(bircol,hetcol,malcol)) + 
  stat_summary(fun.data = function(x) return(c(ymin = mean(x) - 2 * sd(x)/sqrt(length(x)), ymax = mean(x) + 2 * sd(x)/sqrt(length(x)))), geom = "errorbar", color=c(bircol,hetcol,malcol), width=0.15,size=0.8)  + 
  labs(x = "Genotype", y = expression(atop(Maximum~Respiration,(paste(pmol~O[2]~s^{"-1"}~mL^{-1}~mg^{-1}))))) +
  scale_x_discrete(labels = c("*X. birchmanni*", "hybrid", "*X. malinche*"))
maxres
ggsave("output_files/Fig2_max_respiration.pdf",
       maxres, width = 3.75, height = 3, units = "in")


######
# Time to peak post-succinate
######

timetosucc_plot <- ggplot(data, aes(x = species, y = s_time)) +
  theme_bw() +
  ylim(0,900) +
  labs(x = "Genotype", y = "Seconds to Peak CI + CII") +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_markdown()) +
  geom_jitter(width=0.1,col="gray") +
  scale_x_discrete(labels = c("*X. birchmanni*", "hybrid", "*X. malinche*")) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), color=c(bircol,hetcol,malcol), size=0.8)
timetosucc_plot
ggsave("output_files/time_to_succinate_peak_v2.pdf",
       timetosucc_plot, width = 4, height = 4, units = "in")


## Stats analysis
library(car)

norev <- filter(datum, species != "F1-hybrid-reverse") %>%
  mutate(species = as.factor(species))

# Set up contrasts - essentially want to compare hybrids against non-hybrids, while controlling for date
contrast1 <- c(-1/2,1,-1/2)
contrasts(norev$species) = cbind(contrast1)

norevCIlm <- aov(FCF_value ~ species + date, data = filter(norev, FCF == "CI_eff"))
summary.lm(norevCIlm)

norevMRlm <- aov(FCF_value ~ species + date, data = filter(norev, FCF == "Max_prot_norm_resp"))
summary.lm(norevMRlm)

norevETSlm <- aov(FCF_value ~ species + date, data = filter(norev, FCF == "ETS_capacity"))
summary.lm(norevETSlm)

norevCICIIlm <- aov(FCF_value ~ species + date, data = filter(norev, FCF == "CI_CII_eff"))
summary.lm(norevCICIIlm)

norevCIVlm <- aov(FCF_value ~ species + date, data = filter(norev, FCF == "CIV"))
summary.lm(norevCIVlm)

norevLEAKlm <- aov(FCF_value ~ species + date, data = filter(norev, FCF == "Leak_prot_norm_resp"))
summary.lm(norevLEAKlm)

# And now stats on time to peak respiration
data$species <- as.factor(data$species)
contrast1 <- c(-1/2,1,-1/2)
contrasts(data$species) = cbind(contrast1)
# Post-ADP
CItimelm <- aov(adp_time ~ species + date, data = data)
summary.lm(CItimelm)
# Post-succinate
CIItimelm <- aov(s_time ~ species + date, data = data)
summary.lm(CIItimelm)
