library(lubridate) 
library(sjPlot)
library(sjmisc)
library(ggeffects)
library(cowplot)
library(tidyverse)

# load the oxygen consumption data, where the main variable mO2 is oxygen consumption in pmol/min
F2resp <- read.csv("input_files/loligo_F2_embryo_resprate_final_blanksubtracted.csv") %>%
  mutate(Mother = as.factor(Mother),
         Batch = as.factor(interaction(Mother,Batch)),
         Embryo = as.factor(Embryo))#,
         #Date = as.Date(Date,format = "%m/%d/%y"))d

# load the genotype file, which also includes the respirometry well so that the two can be tied together
genotypes1 <- read.csv("input_files/loligo_F2_embryo_batch1_genotypes_dryad.csv") %>%
  mutate(Mother = as.factor(Mother),
         Embryo = as.factor(as.numeric(sapply(sapply(Individual, function(x) strsplit(x, split = "[_-]E")[[1]][2]), function(x) strsplit(x, split = "\\.")[[1]][1]))),
         ndufs5 = as.factor(ndufs5),
         ndufa13 = as.factor(ndufa13),
         chr15 = as.factor(chr15),
         mito = as.factor(mito))
genotypes2 <- read.csv("input_files/loligo_F2_embryo_batch2_genotypes_dryad.csv") %>%
  mutate(Mother = as.factor(Mother),
         Embryo = as.factor(as.numeric(sapply(sapply(Individual, function(x) strsplit(x, split = "[_-]E")[[1]][2]), function(x) strsplit(x, split = "\\.")[[1]][1]))),
         ndufs5 = as.factor(ndufs5),
         ndufa13 = as.factor(ndufa13),
         chr15 = as.factor(chr15),
         mito = as.factor(mito)) %>%
  filter(X100K_reads)
genotypes <- rbind(genotypes1, genotypes2)

#Both genes match HWE in utero:
chisq.test(table(genotypes$ndufa13), p = c(.25, .5, .25))
chisq.test(table(genotypes$ndufa13, genotypes$ndufs5), p = c(.25, .5, .25) %*% t(c(.25, .5, .25)))

# load the morphometrics, calculate the upper quartile of length for each brood (useful for plotting developmental stage later)
# and join to the genotypes based on the mother & embryo numbers
morphometrics <- read.csv("input_files/loligo_F2_embryo_phenotypes.csv") %>%
  mutate(Mother = as.factor(Mother),
         Embryo = as.factor(Embryo),
         Batch = as.factor(interaction(Mother,Batch)),
         yolk_volume = 4/3 * pi * yolk_length_mm * yolk_width_mm * yolk_height_mm / 8,
         yolk_avgdiam = (yolk_length_mm + yolk_width_mm + yolk_height_mm)/3) %>%
  group_by(Mother) %>%
  mutate(brood_upper_quartile = quantile(length_mm, probs = .75, na.rm = T)) %>%
  left_join(genotypes)

morphometrics$ndufs5  <- recode(morphometrics$ndufs5, `0` = "birchmanni", `1` = "heterozygous",`2` = "malinche")
morphometrics$ndufa13  <- recode(morphometrics$ndufa13, `0` = "birchmanni", `1` = "heterozygous",`2` = "malinche")


allresp <- F2resp %>%
  left_join(genotypes)
all_rote <- filter(allresp, Rotenoned == T) %>%
  select(Well, Mother, mO2) %>%
  rename(rotenone_rate = mO2)

# Leaving out brood 3 because they were just fertilized 
# and didn't have much in the way of measurable respiration or interesting phenotypes
F2resp_rate <- filter(allresp, Mother != "3")
nob3_rote <- filter(all_rote, Mother != "3")

F2resp_rate$ndufs5  <- recode(F2resp_rate$ndufs5, `0` = "birchmanni", `1` = "heterozygous",`2` = "malinche")
F2resp_rate$ndufa13  <- recode(F2resp_rate$ndufa13, `0` = "birchmanni", `1` = "heterozygous",`2` = "malinche")

alldata <- left_join(F2resp_rate, morphometrics) %>%
  ungroup()


# This gives the residual effects of all other variables on respiration rate after controlling for length
# (an interesting way to make plots later)
length_controlled <- lm(mO2 ~ length_mm, data = alldata)$residuals
lengthcontrolled_data <- alldata %>%
  filter(!(is.na(length_mm)), !(is.na(mO2))) %>%
  mutate(length_controlled = length_controlled)

# The proportional change in respiration with the application of rotenone
# Turns out to not be that interesting because they all more or less stop breathing with the rotenone
F2resp_diff <- filter(F2resp_rate, Rotenoned == F, !(is.na(ndufs5) & is.na(ndufa13))) %>%
  left_join(nob3_rote) %>%
  mutate(prop_rotenone_sensitive = (mO2 - rotenone_rate) / mO2)

rotemorpho <- left_join(F2resp_diff, morphometrics)
length_controlled_diff <- lm(prop_rotenone_sensitive ~ length_mm, data = filter(rotemorpho,0 < prop_rotenone_sensitive, prop_rotenone_sensitive < 1))$residuals
lengthcontrolled_diffdata <- F2resp_diff %>%
  filter(!(is.na(prop_rotenone_sensitive)), mO2 > 0, 0 < prop_rotenone_sensitive, prop_rotenone_sensitive < 1, rotenone_rate > 0) %>%
  mutate(length_controlled_diff = length_controlled_diff)

length_controlled_heartrate <- lm(heart_rate_bpm ~ length_mm, data = morphometrics)$residuals
lengthcontrolled_heartdata <- morphometrics %>%
  filter(!(is.na(length_mm)), !(is.na(heart_rate_bpm)))
lengthcontrolled_heartdata <- cbind(lengthcontrolled_heartdata, length_controlled_heartrate = length_controlled_heartrate)
#write_csv(lengthcontrolled_heartdata, "Swordtail Dropbox/Schumer_lab_resources/Project_files/Chromosome13_incompatibility/Data/length_controlled_heartrate_morphometrics_genotypes.csv")

length_controlled_snwidth <- lm(sinuatrium_width_mm ~ length_mm, data = morphometrics)$residuals
lengthcontrolled_snwidthdata <- morphometrics %>%
  filter(!(is.na(length_mm)), !(is.na(sinuatrium_width_mm)))
lengthcontrolled_snwidthdata <- cbind(lengthcontrolled_snwidthdata, length_controlled_snwidth = length_controlled_snwidth)
#write_csv(lengthcontrolled_snwidthdata, "Swordtail Dropbox/Schumer_lab_resources/Project_files/Chromosome13_incompatibility/Data/length_controlled_sinuatriumwidth_morphometrics_genotypes.csv")

length_controlled_headwidth <- lm(head_width_mm ~ length_mm, data = morphometrics)$residuals
lengthcontrolled_headwidthdata <- morphometrics %>%
  filter(!(is.na(length_mm)), !(is.na(head_width_mm)))
lengthcontrolled_headwidthdata <- cbind(lengthcontrolled_headwidthdata, length_controlled_headwidth = length_controlled_headwidth)

length_controlled_yolkdiam <- lm(yolk_avgdiam ~ length_mm, data = filter(morphometrics, Damaged == F))$residuals
lengthcontrolled_yolkdiamdata <- morphometrics %>%
  filter(!(is.na(length_mm)), !(is.na(yolk_avgdiam)), Damaged == F)
lengthcontrolled_yolkdiamdata <- cbind(lengthcontrolled_yolkdiamdata, length_controlled_yolkdiam = length_controlled_yolkdiam)

broodcontrolled_length <- lm(length_mm ~ brood_upper_quartile, data = filter(morphometrics, Damaged == F))$residuals
broodcontrolled_lengthdata <- morphometrics %>%
  filter(!(is.na(length_mm)), !(is.na(brood_upper_quartile)), Damaged == F)
broodcontrolled_lengthdata <- cbind(broodcontrolled_lengthdata, broodcontrolled_length = broodcontrolled_length)

# testing effects of our nuclear genes of interest, ndufs5 and ndufa13, the effect of length, and all their interactions with each other
# on various phenotypes, including Mother as a non-interacting blocking variable to control for maternal/batch effects
### NOTE that we use the morphometrics dataset here, rather than alldata, because alldata includes two rows per individual (pre- and post-rotenone)
### So using alldata without filtering based on rotenone status would mean artificially doubling the morphometrics sample sizes
summary(aov(length_mm ~ ndufs5 * ndufa13 * brood_upper_quartile + Mother, data = filter(morphometrics, Damaged == F)))
summary(aov(heart_rate_bpm ~ ndufs5 * ndufa13 * length_mm + Mother, data = filter(morphometrics, Damaged == F)))
summary(aov(sinuatrium_width_mm ~ ndufs5 * ndufa13 * length_mm + Mother, data = filter(morphometrics, Damaged == F)))
summary(aov(head_width_mm ~  ndufs5 * ndufa13 * length_mm + Mother, data = filter(morphometrics, Damaged == F)))
summary(aov(yolk_avgdiam ~ ndufs5 * ndufa13 * length_mm + Mother, data = filter(morphometrics, Damaged == F)))


# Now a model to test the effect of the same variables on oxygen consumption
# this time using filtered alldata as the input data because we need the respiration columns
# and batch as the blocking variable, since in one case I had to split the brood across two respirometry plates,
# and we want to be sure to control for differences between plates, not just between Mothers
summary(ratemodel) <- aov(mO2 ~ ndufs5 * ndufa13 * length_mm + Batch, data = filter(alldata, Rotenoned == F, Damaged == F))

summary(aov(prop_rotenone_sensitive ~ ndufs5 * ndufa13 * length_mm + Batch, data = filter(rotemorpho, 0 < prop_rotenone_sensitive, prop_rotenone_sensitive < 1, Damaged == F)))
TukeyHSD(aov(length_controlled_diff ~ ndufs5, data = filter(lengthcontrolled_diffdata, 0 < prop_rotenone_sensitive, prop_rotenone_sensitive < 1)))



## Plotting
malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

theme_set(theme_bw())
ratemodel2 <- aov(mO2 ~ ndufs5 * ndufa13 * length_mm + Batch, data = filter(alldata, Rotenoned == F, ndufs5 != "birchmanni"))
df <- ggpredict(ratemodel2, terms = c("length_mm","ndufs5","ndufa13"))
marginal <- ggplot(df, aes(x = x, y = predicted)) + 
  ggtitle("ndufa13") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~facet) +
  geom_line(aes(color=group)) +
  scale_colour_manual(values = c(hetcol, malcol, bircol)) +
  scale_fill_manual(values = c(hetcol, malcol, bircol)) + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.4) +
  labs(x = "Length (cm)", y = expression(atop(Predicted~O[2]~"Consumption Rate", "(pmol / min, Mean ± 2 SE)"))) +
  guides(fill=guide_legend(title="ndufs5"), colour=guide_legend(title="ndufs5")) +
  coord_cartesian(x = c(2.5, 11))
  marginal
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Chromosome13_incompatibility/Figures/ndufs5_ndufa13_length_o2_marginal_plot.pdf",
              marginal, width = 8, height = 4, units = "in")

# Plot that shows divergence in size based on genotype as development progresses
length_genotype_ndufs5_plot <- ggplot(filter(morphometrics, !is.na(ndufs5), !is.na(ndufa13)), aes(x = brood_upper_quartile, y = length_mm, color = ndufs5, fill = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.25, 0.80),
        legend.box = "horizontal",
        panel.grid = element_blank()) +
  labs(x = "Brood Length (mm)", y = "Individual Length (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm")
length_genotype_ndufs5_plot
length_genotype_ndufa13_plot <- ggplot(filter(morphometrics, !is.na(ndufs5), !is.na(ndufa13)), aes(x = brood_upper_quartile, y = length_mm, color = ndufa13, fill = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.25, 0.80),
        legend.box = "horizontal",
        panel.grid = element_blank()) +
  labs(x = "Brood Length (mm)", y = "Individual Length (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm")
length_genotype_ndufa13_plot
length_genotype_twopanel <- cowplot::plot_grid(length_genotype_ndufs5_plot, length_genotype_ndufa13_plot, nrow = 2, labels = "AUTO")
length_genotype_twopanel
ggsave("output_files/stageproxy_ndufs5_ndufa13_length_plot.pdf",
       length_genotype_twopanel, width = 6, height = 8, units = "in")

ndufs5_lengthinbrood_plot <- ggplot(filter(broodcontrolled_lengthdata, !(is.na(ndufs5))), aes(x = ndufs5, y = broodcontrolled_length, fill = ndufs5, color = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.19, 0.80),
        panel.grid = element_blank()) +
  #  coord_cartesian(y = c(-20, 30))+
  labs(x = "ndufs5", y = "Brood-Controlled Length (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_x_discrete(labels = c("X. bir", "Heterozygote", "X. mal")) +
  guides(color = F, fill = F) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_jitter(width = 0.25, color = "gray", alpha = 0.5) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), size=0.8)
ndufs5_lengthinbrood_plot
ggsave("output_files/stagecorrected_ndufs5_length_plot.pdf",
       ndufs5_lengthinbrood_plot, width = 3, height = 3, units = "in")

ndufa13_lengthinbrood_plot <- ggplot(filter(broodcontrolled_lengthdata, !(is.na(ndufa13))), aes(x = ndufa13, y = broodcontrolled_length, fill = ndufa13, color = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.19, 0.80),
        panel.grid = element_blank()) +
#  coord_cartesian(y = c(-20, 30))+
  labs(x = "ndufa13", y = "Brood-Controlled Length (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_x_discrete(labels = c("X. bir", "Heterozygote", "X. mal")) +
  guides(color = F, fill = F) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_jitter(width = 0.25, color = "gray", alpha = 0.5) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), size=0.8)
ndufa13_lengthinbrood_plot
ggsave("output_files/stagecorrected_ndufa13_length_plot.pdf",
       ndufa13_lengthinbrood_plot, width = 3, height = 3, units = "in")

# new heartrate plot controlling for length, so that it can be plotted like a bar chart
# bar-chart style heart rate
heartrate_plot <- ggplot(filter(lengthcontrolled_heartdata, !(is.na(ndufa13))), aes(x = ndufa13, y = length_controlled_heartrate, fill = ndufa13, color = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.19, 0.80),
        panel.grid = element_blank()) +
  coord_cartesian(y = c(-20, 30))+
  labs(x = "ndufa13", y = "Length-Controlled Heart Rate (bpm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_x_discrete(labels = c("X. bir", "Heterozygote", "X. mal")) +
  guides(color = F, fill = F) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_jitter(width = 0.25, color = "gray", alpha = 0.5) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), size=0.8)
heartrate_plot
ggsave("output_files/loligo_ndufa13_heartrate_lengthresiduals_plot.pdf",
       heartrate_plot, width = 3, height = 3, units = "in")
TukeyHSD(aov(length_controlled_heartrate ~ ndufa13, filter(lengthcontrolled_heartdata, !(is.na(ndufa13)))))

# Plot of heartrate with length on the x
heartrate_ndufs5_plot <- ggplot(filter(morphometrics, !(is.na(ndufs5)), !(is.na(ndufa13)), Damaged == F), aes(x = length_mm, y = heart_rate_bpm, fill = ndufs5, color = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.25, 0.80),
        legend.box = "horizontal",
        panel.grid = element_blank()) +
  labs(x = "Individual Length (mm)", y = "Heart Rate (bpm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm")
heartrate_ndufs5_plot
heartrate_ndufa13_plot <- ggplot(filter(morphometrics, !(is.na(ndufs5)), !(is.na(ndufa13)), Damaged == F), aes(x = length_mm, y = heart_rate_bpm, fill = ndufa13, color = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.25, 0.80),
        legend.box = "horizontal",
        panel.grid = element_blank()) +
  labs(x = "Individual Length (mm)", y = "Heart Rate (bpm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm")
heartrate_ndufa13_plot
heartrate_length_twopanel <- cowplot::plot_grid(heartrate_ndufs5_plot, heartrate_ndufa13_plot, nrow = 2, labels = "AUTO")
heartrate_length_twopanel
ggsave("output_files/loligo_ndufs5_ndufa13_heartrate_plot.pdf",
       heartrate_length_twopanel, width = 6, height = 8, units = "in")

# Plot of sinuatrium width with length on the x
snwidth_length_ndufs5_plot <- ggplot(filter(morphometrics, !(is.na(ndufs5))), aes(x = length_mm, y = sinuatrium_width_mm, fill = ndufs5, color = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.80),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Sinu-atrium Width (mm)") +
  geom_point() +
  geom_smooth(method = "lm")
snwidth_length_ndufs5_plot
snwidth_length_ndufa13_plot <- ggplot(filter(morphometrics, !(is.na(ndufa13))), aes(x = length_mm, y = sinuatrium_width_mm, fill = ndufa13, color = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.80),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Sinu-atrium Width (mm)") +
  geom_point() +
  geom_smooth(method = "lm")
snwidth_length_ndufa13_plot
snwidth_length_twopanel <- cowplot::plot_grid(snwidth_length_ndufs5_plot, snwidth_length_ndufa13_plot, nrow = 2, labels = "AUTO")
snwidth_length_twopanel
ggsave("output_files/snwidth_ndufs5_ndufa13_length_plot.pdf",
       snwidth_length_twopanel, width = 6, height = 8, units = "in")

# bar-chart style sinu-atrium width plot
snwidth_plot <- ggplot(filter(morphometrics, !(is.na(ndufa13))), aes(x = ndufa13, y = sinuatrium_width_mm, color = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  guides(color = F) +
  labs(x = "ndufa13", y = "Sinu-atrium Width (mm)") +
  geom_jitter(width = 0.25, color = "gray") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), size=0.8)
snwidth_plot
ggsave("output_files/loligo_ndufa13_sinuatrium_width_plot.pdf", 
       snwidth_plot, width = 4, height = 4, units = "in")
TukeyHSD(aov(sinuatrium_width_mm ~ ndufa13, filter(morphometrics, !(is.na(ndufa13)))))

snwidth_plot2 <- ggplot(filter(lengthcontrolled_snwidthdata, !(is.na(ndufa13))), aes(x = ndufa13, y = length_controlled_snwidth, color = ndufa13)) +
  theme_bw() +
  coord_cartesian(y = c(-.1, .1)) +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_x_discrete(labels = c("X. bir", "Heterozygote", "X. mal")) +
  guides(color = F) +
  labs(x = "ndufa13", y = expression(atop("Length-Controlled", "Sinu-atrium Width (mm)"))) +
  geom_jitter(width = 0.25, color = "gray", alpha = 0.5) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), size=0.8)
snwidth_plot2
ggsave("output_files/loligo_ndufa13_sinuatrium_width_lenthcontrolled_plot.pdf",
       snwidth_plot2, width = 3, height = 3, units = "in")
TukeyHSD(aov(length_controlled_snwidth ~ ndufa13, filter(lengthcontrolled_snwidthdata, !(is.na(ndufa13)))))

# head width, length, and ndufs5 + ndufa13 genotype
headwidth_ndufs5_length_plot <- ggplot(filter(morphometrics, !is.na(ndufs5)), aes(x = length_mm, y = head_width_mm, color = ndufs5, fill = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Head Width (mm)") +
  geom_point() +
  geom_smooth(method = "lm")
headwidth_ndufs5_length_plot
headwidth_ndufa13_length_plot <- ggplot(filter(morphometrics, !is.na(ndufa13)), aes(x = length_mm, y = head_width_mm, color = ndufa13, fill = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Head Width (mm)") +
  geom_point() +
  geom_smooth(method = "lm")
headwidth_ndufa13_length_plot
headwidth_genotype_twopanel <- cowplot::plot_grid(headwidth_ndufs5_length_plot, headwidth_ndufa13_length_plot, nrow = 2, labels = "AUTO")
headwidth_genotype_twopanel
ggsave("output_files/headwidth_ndufs5_ndufa13_length_plot.pdf",
       headwidth_genotype_twopanel, width = 6, height = 8, units = "in")

lengthcontrolled_headwidth_ndufs5_plot <- ggplot(filter(lengthcontrolled_headwidthdata, !is.na(ndufs5)), aes(x = ndufs5, y = length_controlled_headwidth, color = ndufs5)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "ndufs5", y = "Length-Controlled Head Width (mm)") +
  guides(color = F) +
  geom_jitter(width = 0.2, color = "gray") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))
lengthcontrolled_headwidth_ndufs5_plot
ggsave("output_files/ndufs5_head_width_lenthcontrolled_plot.pdf",
       lenthcontrolled_headwidth_ndufs5_plot, width = 6, height = 5, units = "in")


yolk_stage_ndufs5_plot <- ggplot(filter(morphometrics, !is.na(ndufs5), Damaged == F), aes(x = brood_upper_quartile, y = yolk_avgdiam, color = ndufs5, fill = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.80),
        panel.grid = element_blank()) +
  labs(x = "Brood Length (mm)", y = "Yolk Diameter (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm") 
yolk_stage_ndufs5_plot
yolk_stage_ndufa13_plot <- ggplot(filter(morphometrics, !is.na(ndufa13), Damaged == F), aes(x = brood_upper_quartile, y = yolk_avgdiam, color = ndufa13, fill = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.80),
        panel.grid = element_blank()) +
  labs(x = "Brood Length (mm)", y = "Yolk Diameter (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm") 
yolk_stage_ndufa13_plot
yolk_genotype_twopanel <- cowplot::plot_grid(yolk_stage_ndufs5_plot, yolk_stage_ndufa13_plot, nrow = 2, labels = "AUTO")
yolk_genotype_twopanel
ggsave("output_files/ndufs5_ndufa13_stageproxy_yolkdiam_plot.pdf",
       yolk_genotype_twopanel, width = 6, height = 8, units = "in")


# Yolk diameter as a function of length and genotype
yolk_length_ndufs5_plot <- ggplot(filter(morphometrics, !is.na(ndufs5), Damaged == F), aes(x = length_mm, y = yolk_avgdiam, color = ndufs5, fill = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.80),
        panel.grid = element_blank()) +
  labs(x = "Individual Length (mm)", y = "Yolk Diameter (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm")
yolk_length_ndufs5_plot
yolk_length_ndufa13_plot <- ggplot(filter(morphometrics, !is.na(ndufa13), Damaged == F), aes(x = length_mm, y = yolk_avgdiam, color = ndufa13, fill = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.80),
        panel.grid = element_blank()) +
  labs(x = "Individual Length (mm)", y = "Yolk Diameter (mm)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  geom_point() +
  geom_smooth(method = "lm")
yolk_length_ndufa13_plot
yolk_length_twopanel <- cowplot::plot_grid(yolk_length_ndufs5_plot, yolk_length_ndufa13_plot, nrow = 2, labels = "AUTO")
ggsave("output_files/ndufs5_ndufa13_length_yolkdiam_plot.pdf",
       yolk_length_twopanel, width = 6, height = 8, units = "in")

lengthcontrolled_yolk_ndufs5plot <- ggplot(filter(lengthcontrolled_yolkdiamdata, !is.na(ndufs5)), aes(x = ndufs5, y = length_controlled_yolkdiam)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  #  facet_wrap(~ Rotenoned) +
  labs(x = "ndufs5 Genotype", y = "Length-Controlled Yolk Diameter (mm)")  +
  geom_jitter(width=0.1,col="gray") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), color=c( bircol,hetcol, malcol), size=0.8)
lengthcontrolled_yolk_ndufs5plot
ggsave("output_files/ndufs5_yolk_lenthcontrolled_plot.pdf",
       lengthcontrolled_yolk_ndufs5plot, width = 6, height = 5, units = "in")
TukeyHSD(aov(length_controlled_yolkdiam ~ ndufs5, data = filter(lengthcontrolled_yolkdiamdata)))



# Beautifying the respirometry data for plotting
plottable <- mutate(F2resp_rate, Rotenoned = ifelse(Rotenoned, "Post-Rotenone", "Pre-Rotenone"))
plottable$Rotenoned <- factor(plottable$Rotenoned, levels = c("Pre-Rotenone", "Post-Rotenone"))

# Plot just the effects of ndufs5 on mO2
ndufs5plot = ggplot(filter(plottable, !is.na(ndufs5), !is.na(Rotenoned)), aes(x = ndufs5, y = mO2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_jitter(width=0.1,col="gray") +
  labs(x = "ndufs5 Genotype", y = "Oxygen Consumption Rate\n(pmol / min, Mean ± 2 SE)") +
  scale_x_discrete(labels = c("X. birchmanni", "heterozygous", "X. malinche")) +
  facet_grid(~Rotenoned) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), color=c(bircol,hetcol, malcol, bircol,hetcol,malcol), size=0.8)
ndufs5plot
ggsave("output_files/loligo_ndufs5_plot.pdf",
      ndufs5plot, width = 6, height = 5, units = "in")  

# Plot just the effects of ndufs5 on mO2 pre-rotenone
ndufs5plot_preonly = ggplot(filter(plottable, !is.na(ndufs5), Rotenoned == "Pre-Rotenone"), aes(x = ndufs5, y = mO2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_jitter(width=0.1,col="gray") +
  labs(x = "ndufs5 Genotype", y = "Oxygen Consumption Rate\n(pmol / min, Mean ± 2 SE)") +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  #facet_grid(~Rotenoned) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), color=c(bircol,hetcol, malcol), size=0.8)
ndufs5plot_preonly
ggsave("output_files/loligo_ndufs5_prerotenone_plot.pdf",
       ndufs5plot_preonly, width = 3.75, height = 3, units = "in")  

# oxygen consumption, length, and ndufs5 genotype
o2_ndufs5_length_plot <- ggplot(filter(alldata, Rotenoned == F, !(is.na(ndufs5))), aes(x = length_mm, y = mO2, color = ndufs5, fill = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.19, 0.80),
        panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Oxygen Consumption Rate\n(pmol / min, Mean ± 2 SE)") +
  geom_point() +
  geom_smooth(method = "lm")
o2_ndufs5_length_plot
ggsave("output_files/loligo_ndufs5_o2_length_plot.pdf",
       o2_ndufs5_length_plot, width = 6, height = 5, units = "in")

# oxygen consumption, length, and ndufs5 + ndufa13 genotype
o2_ndufs5_length_plot <- ggplot(filter(alldata, !is.na(ndufs5), Rotenoned == F), aes(x = length_mm, y = mO2, color = ndufs5, fill = ndufs5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.5, 0.8)) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Oxygen Consumption Rate\n(pmol / min, Mean ± 2 SE)") +
  geom_point() +
  geom_smooth(method = "lm")
o2_ndufs5_length_plot
o2_ndufa13_length_plot <- ggplot(filter(alldata, !is.na(ndufa13), Rotenoned == F), aes(x = length_mm, y = mO2, color = ndufa13, fill = ndufa13)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.5, 0.8)) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "Individual Length (mm)", y = "Oxygen Consumption Rate\n(pmol / min, Mean ± 2 SE)") +
  geom_point() +
  geom_smooth(method = "lm")
o2_ndufa13_length_plot
o2_length_twopanel <- cowplot::plot_grid(o2_ndufs5_length_plot, o2_ndufa13_length_plot, nrow = 2, labels = "AUTO")
o2_length_twopanel
ggsave("output_files/loligo_ndufs5_ndufa13_o2_length_plot.pdf",
       o2_length_twopanel, width = 6, height = 8, units = "in")

# oxygen consumption and ndufs5 + ndufa13 genotype (not controlling for length)
o2_ndufs5_ndufa13_plot <- ggplot(filter(alldata, !is.na(ndufs5), !is.na(ndufa13), Rotenoned == F), aes(x = interaction(ndufa13,ndufs5), y = mO2, color = ndufa13)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "ndufs5", y = expression(atop(O[2]~"Consumption Rate", "(pmol / min, Mean ± 2 SE)"))) +
  geom_jitter(width = 0.2, color = "gray") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))
o2_ndufs5_ndufa13_plot
ggsave("output_files/loligo_ndufs5_ndufa13_o2_plot.pdf",
       o2_ndufs5_ndufa13_plot, width = 6, height = 5, units = "in")

# oxygen consumption and ndufs5 + ndufa13 genotype (controlling for length)
o2_ndufs5_ndufa13_plot <- ggplot(filter(lengthcontrolled_data, !is.na(ndufs5), !is.na(ndufa13), Rotenoned == F), aes(x = interaction(ndufa13,ndufs5), y = length_controlled, color = ndufa13)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c(bircol,hetcol, malcol)) +
  labs(x = "ndufs5", expression(atop("Length-Controlled"~O[2]~"Consumption", "(pmol / min, Mean ± 2 SE)"))) +
  geom_jitter(width = 0.2, color = "gray") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))
o2_ndufs5_ndufa13_plot


# oxygen consumption and ndufs5 genotype (controlling for length)
lengthcontrolled_o2_ndufs5plot <- ggplot(filter(lengthcontrolled_data, Rotenoned == F, !is.na(ndufs5)), aes(x = ndufs5, y = length_controlled)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #coord_cartesian(y = c(-100, 275)) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  guides(fill = F) +
  labs(x = "ndufs5", y = expression(atop("Length-Controlled"~O[2]~"Consumption", "(pmol / min, Mean ± 2 SE)"))) +
  scale_x_discrete(labels = c("X. bir", "Heterozygote", "X. mal")) +
  geom_jitter(width=0.1,col="gray", alpha = 0.5) +
  #geom_violin(aes(fill = ndufs5), alpha = 0.25)+
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), color=c( bircol,hetcol, malcol), size=0.55)
lengthcontrolled_o2_ndufs5plot
ggsave("output_files/loligo_ndufs5_mO2_lengthresiduals_plot.pdf",
       lengthcontrolled_o2_ndufs5plot, width = 3, height = 3, units = "in")
TukeyHSD(aov(length_controlled ~ ndufs5, data = filter(lengthcontrolled_data, Rotenoned == F)))

# Proportion rotenone-sensitive respiration plot
propsensitive_ndufs5plot <- ggplot(filter(F2resp_diff, 0 < prop_rotenone_sensitive, prop_rotenone_sensitive < 1), aes(x = ndufs5, y = prop_rotenone_sensitive)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c(bircol,hetcol, malcol)) +
  guides(fill = F) +
  labs(x = "ndufs5", y = expression(atop("Proportional Respiration Decline Post-Rotenone", "(pmol / min, Mean ± 2 SE)"))) +
  scale_x_discrete(labels = c("X. bir", "Heterozygote", "X. mal")) +
  geom_jitter(width=0.1,col="gray", alpha = 0.5) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)), color=c( bircol,hetcol, malcol), size=0.55)
propsensitive_ndufs5plot
ggsave("output_files/loligo_ndufs5_propsensitive_plot.pdf",
       propsensitive_ndufs5plot, width = 6, height = 5, units = "in")
TukeyHSD(aov(length_controlled_diff ~ ndufs5, data = filter(lengthcontrolled_diffdata, Rotenoned == F)))
