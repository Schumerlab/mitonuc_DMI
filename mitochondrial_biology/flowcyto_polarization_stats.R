library(tidyverse)

fc_data <- read.csv('Data_Flowcytometry_xipho.csv') %>%
  mutate(diff = Untreated - FCCP,
         prop_diff = FCCP/Untreated)
fc_data$Genotype <- as.factor(fc_data$Genotype) 

fc_data_long <- pivot_longer(fc_data, Untreated:FCCP, names_to = "Treatment", values_to = "fluorescence")
fc_data_long$Genotype <- as.factor(fc_data_long$Genotype)

datasumm <- group_by(fc_data_long, Treatment, Genotype) %>%
  summarize(mean = mean(fluorescence),
            stdev = sd(fluorescence),
            se = stdev/sqrt(n()))

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

fc_plot <- ggplot(datasumm, aes(x = Genotype, y = mean, fill = Treatment, color = Treatment)) +
  theme_minimal() +
#  theme(legend.position = "none") +
  geom_point(stat = "identity", cex = 2.5) +
  geom_jitter(data = fc_data_long, aes(y = fluorescence),  width=0.1, alpha = 0.5) +
  labs(x = "Genotype", y = "Median Fluorescence Intensity") +
  scale_color_manual(values = c(hetcol, rgb(50/255,50/255,50/255))) +
  scale_y_log10() +
  geom_errorbar(aes(ymin = mean - 2 * se, ymax = mean + 2*se), width=0.15,size=0.8)
fc_plot
ggsave('~/Library/CloudStorage/Box-Box/Schumer_lab_resources/Project_files/Chromosome13_incompatibility/Figures/flowcytometry_graph_log.pdf',
       fc_plot, width = 6, height = 4, units = "in")

summary(aov(log(Untreated) ~ Genotype + FCCP, data = fc_data))

fc_data <- fc_data %>%
  mutate(diff = Untreated - FCCP,
         prop_diff = FCCP/Untreated,
         arcsin_prop_diff = asin(sqrt(prop_diff)))

summary(aov(1 - prop_diff ~ Genotype, data = fc_data))

