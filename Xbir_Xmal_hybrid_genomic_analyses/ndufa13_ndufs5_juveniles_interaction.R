library(tidyverse)
library(ggmosaic)
juvenilesMar2023 <- read.csv("input_files/genotypes_morphometrics_F2_Mar2023.csv") %>%
  mutate(Individual = id,
         ndufa13 = recode(ndufa13, "0" = "birchmanni", "1" = "heterozygous", "2" = "malinche"),
         ndufs5 = recode(ndufs5, "0" = "birchmanni", "1" = "heterozygous", "2" = "malinche")) %>% 
  select(Individual, ndufs5, ndufa13, chr15, mito)
juvie1 <- table(select(juvenilesMar2023,ndufa13,ndufs5))
juvenilesJul2023 <- read.csv("Documents/schumer_lab/mitonuclear_incompatibility/juvenile_F2_Oct22_genotypes.csv") %>%
  mutate(ndufa13 = recode(ndufa13, "0" = "birchmanni", "1" = "heterozygous", "2" = "malinche"),
         ndufs5 = recode(ndufs5, "0" = "birchmanni", "1" = "heterozygous", "2" = "malinche")) 
juvie2 <- table(select(juvenilesJul2023,ndufa13,ndufs5))

F2_names <- read_tsv("input_files/genotypes_outputfile_allchrs_allF2_orlater_stocktanks_March2019_removeF1s.tsv_chr6_chr13_withIDs")
F2s <- read_tsv("input_files/combined_focal_region_genotypes_outputfile_allchrs_allF2_orlater_stocktanks_March2019_removeF1s.tsv_chr6_chr13") %>%
  mutate(Individual = F2_names$id, 
         ndufa13 = recode(`ScyDAA6-2393-HRSCAF-2888:12131248`, "0" = "birchmanni", "1" = "heterozygous", "2" = "malinche"),
         ndufs5 = recode(`ScyDAA6-1934-HRSCAF-2318:2066737`, "0" = "birchmanni", "1" = "heterozygous", "2" = "malinche"),
         chr15 = NA, mito = 2) %>%
  select(Individual, ndufs5, ndufa13, chr15, mito)

allfish <- rbind(F2s, juvenilesMar2023, juvenilesJul2023)

table(F2_names$`ScyDAA6-2393-HRSCAF-2888:12131248`)

juvies = juvie1+juvie2
apply(juvies, 1, sum)
7/sum(7,41,26)


juvie_mean = 7/74
juvie_se = sqrt(7/74 * (1- 7/74)/74)
embryo_mean = 46/208
embryo_se = sqrt(46/208 * (1- 46/208)/208)
adult_mean = 28/932
adult_se = sqrt(28/932 * (1- 28/932)/932)


freq_over_age <- data.frame(age = c("embryo", "juvenile", "adult"), 
                            freq_est = c(embryo_mean, juvie_mean, adult_mean),
                            freq_se = c(embryo_se, juvie_se, adult_se)) %>%
  mutate(age = factor(age, levels = c("embryo","juvenile","adult")))

freq_age_plot <- ggplot(freq_over_age, aes(x = age, y = freq_est)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_point(color = bircol) +
  labs(x = "Age", y = expression(atop("Frequency Homozygous X. birchmanni", "(Estimate ± 2 SE)"))) +
  scale_x_discrete(labels = c("Embryo", expression(atop("Juvenile","(3-5 mo)")), expression(atop("Adult","(6+ mo)")))) +
  geom_pointrange(aes(ymin = freq_est - 2 * freq_se, ymax = freq_est + 2 * freq_se), color = bircol)
freq_age_plot
ggsave("output_files/ndufa13_frequency_age_barplot.pdf",
       freq_age_plot, width = 3, height = 3, units = "in")

#Test for deviation from HWE in juveniles
chisq.test(c(7,41, 26), p = c(.25, .5, .25))
# Frequency is 9.46%, compared to 22.12% in utero

F2_names %>%
  select(`ScyDAA6-2393-HRSCAF-2888:12131248`) %>%
  drop_na() %>%
  group_by(`ScyDAA6-2393-HRSCAF-2888:12131248`) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# Frequency of birchmanni ndufa13 in adults is 3.00%
# SE should be 

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

adult <- table(select(F2s,ndufa13, ndufs5))[1:3, 2:3]
apply(adult, 1, sum)[1]/sum(apply(adult1, 1, sum))
final <- juvies + adult 
chisq.test(final)
chisq.test(juvies, simulate.p.value = T, B = 1000000)

neutral <- c(1/3, 2/3)
selection <- c((1 * .09), (2 * (.09) * (1 - .996 * .12))) / sum((1 * .09), (2 * (.09) * (1 - .996 * .12))) 
# ^ ndufa13 homozygote vs. heterozygote freq * ndufs5 selection's effect on a homozygotes vs. heterozygotes
expected_allndufa13 <- rev(apply(final, 2, sum))/sum(final)
observed <- c(18/33, 15/33)


ratios <- rbind(neutral, selection, expected_allndufa13, observed) 
ratios <- data.frame(ratios, rownames(ratios))
colnames(ratios) <- c("homo", "het", "scenario")
ratios <- ratios %>%
  pivot_longer(homo:het, names_to = "genotype", values_to = "count")
ratios$scenario <- factor(ratios$scenario, levels = c("neutral", "selection", "expected_allndufa13", "observed"))

independent <- as.table(rbind(c((1 - .996) * (1 - .91), 2 * (1 - .996 * .12) * (1 - .91), (1 - .91)), 
                              c(2 * (1 - .996) * (1 - .91 * .09), 4 * (1 - .996 * .12) * (1 - .91 * .09), 2 * (1 - .91 * .09)),
                              c(1 - .996, 2 * (1 - .996 * .12), 1)))

names(dimnames(independent)) = c("ndufa13", "ndufs5")
rownames(independent) = c("bir","het","mal")
colnames(independent) = c("bir","het","mal")
independent
no_ndufs5incomp_expect <- independent[,2:3]
# Test if ndufs5 and ndufa13 are independent
chisq.test(final)

# Getting likelihood of observed ndufs5 counts from the empirical distribution:
set.seed(15)
sampler <- filter(allfish, ndufs5 != "birchmanni")
emp_ratios <- lapply(1:10000, function(x) table(rev(sample(sampler$ndufs5, 33, replace = F)))) %>%
  bind_rows()
tester <- ecdf(emp_ratios$heterozygous)
tester(15)


# Test if ndufs5 and ndufa13 counts differ from those expected for each genotype in case of independent selection
# (MES Please add)



## Plotting:

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

colorscale <- c("homo" = malcol, "het" = hetcol)

plot_props <- data.frame(ndufa13 = row.names(final), 
                         hetcount = final[,1], 
                         malcount = final[,2], 
                         total = apply(final, 1, sum)) %>%
  mutate(hetfreq = hetcount/total,
         hetfreq_se = sqrt(hetfreq * (1 - hetfreq)/total))
                         
hetsonly_plot <- ggplot(plot_props, aes(x = ndufa13, y = hetfreq, color = ndufa13)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_markdown()) +
  coord_cartesian(y = c(0,1)) +
  scale_color_manual(values = c(bircol, hetcol, malcol)) +
  labs(x = expression(italic("ndufa13")), y = expression(atop("Frequency", paste("Heterozygous ", italic("ndufs5"))))) +
  guides(color = F) +
  geom_pointrange(aes(ymin = hetfreq - 2 * hetfreq_se, ymax = hetfreq + 2 * hetfreq_se), size = 0.8) +
  geom_hline(yintercept = 2 * (1 - .996 * .12) * (1 - .91)/ sum(2 * (1 - .996 * .12) * (1 - .91), (1 - .91)), lty = "dashed") +
  geom_text(x = c(1,2,3), y = c(0.1,0.1,0.1), label = c("N = 33", "626", "341")) +
  scale_x_discrete(labels = c("*X. bir*", "Heterozygote", "*X. mal*")) +
  geom_point()
hetsonly_plot
ggsave("output_files/ndufs5het_depletion_plot.pdf",
       hetsonly_plot, width = 3, height = 3, units = "in")