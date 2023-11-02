library(tidyverse)

data_pep <- read.csv("input_files/PRM_percentmalinche_WLLPQSGEQPYK.csv",
                 header = T) %>%
  filter(intervals == "original") %>%
  pivot_longer(ind1:ind5, names_to = "indiv", values_to = "percent_malinche") %>%
  mutate(isotope = recode(peptide, WLLPQSGEQPYK = "endogenous", WLLPQSGEQPYK_heavy = "heavy_raw", WLLPQSGEQPYK_heavy_corrected = "heavy"))

t.test(filter(data_pep, isotope == "endogenous")$percent_malinche, mu = 0.5)
t.test(filter(data_pep, isotope == "heavy")$percent_malinche, mu = 0.5)

datasumm <- group_by(data_pep, isotope) %>%
  summarize(mean = mean(percent_malinche),
            stdev = sd(percent_malinche),
            se = stdev/sqrt(n()))

peptide_imbalance_plot <- ggplot(filter(datasumm, isotope %in% c("endogenous", "heavy")), aes(x = isotope, y = mean, fill = isotope, color = isotope)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(panel.grid = element_blank()) +
  geom_point(stat = "identity") +
  geom_jitter(data = filter(data_pep, isotope %in% c("endogenous", "heavy")), aes(y = percent_malinche), width = 0.2, color = "grey") +
  lims(y = c(0.35, 0.65)) +
  labs(x = "NDUFS5 Peptide", y = expression(paste("Proportion ",italic("X. malinche")))) +
  scale_color_manual(values = c(hetcol, rgb(50/255,50/255,50/255), rgb(150/255,150/255,10/255))) +
  scale_x_discrete(labels = c("Endogenous", "Spike-In")) +
  geom_hline(aes(yintercept = 0.5), lty = 2, color = "grey") + 
  geom_errorbar(aes(ymin = mean - 2 * se, ymax = mean + 2*se), width=0.15,size=0.8)
peptide_imbalance_plot
