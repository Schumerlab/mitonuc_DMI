library(tidyverse)
library("ppcor")
theme_set(theme_bw())
chr2length <- 31064592

mother_embryo_stage_genotypes <-read.csv(file="ndufs5_ndufa13_chr15_mito_stages.csv",sep=",",head=TRUE)

# Test of effects of ndufs5 genotype on developmental lag:
t.test(filter(mother_embryo_stage_genotypes, mitotype == "malinche", ndufs5 != "birchmanni", stage > 3.5)$relative_stage,
       filter(mother_embryo_stage_genotypes, mitotype == "malinche", ndufs5 == "birchmanni", stage > 3.5)$relative_stage)

t.test(filter(mother_embryo_stage_genotypes, mitotype == "malinche", ndufa13 != "birchmanni", stage > 3.5)$relative_stage,
       filter(mother_embryo_stage_genotypes, mitotype == "malinche", ndufa13 == "birchmanni", stage > 3.5)$relative_stage)

################
###embryo stage vs. ndufs5 genotype figure
###############


malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

binned_chr13_plot <- ggplot(filter(mother_embryo_stage_genotypes, !is.na(ndufs5), stage < 12, mitotype == "malinche"), aes(x = stage, fill = ndufs5, group = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.70),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c(bircol, hetcol, malcol)) +
  guides(fill=guide_legend(title=expression(paste(italic("ndufs5")," Genotype"))))+
  labs(x = "Stage", y = "Count") +
  geom_histogram(binwidth = 1.5)
binned_chr13_plot
ggsave("output_files/CALL_malmito_chr13_embryos_hist.pdf", binned_chr13_plot, 
       device = "pdf", width = 5, height = 2.92)


################
###embryo relative stage vs. genotype figure
###############

binned_chr13_relative <- ggplot(filter(mother_embryo_stage_genotypes, !is.na(ndufs5), stage < 12, mitotype == "malinche"), aes(x = relative_stage, fill = ndufs5, group = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.70),
        panel.grid = element_blank()) +
  #scale_x_reverse() +
  scale_fill_manual(values = c(bircol, hetcol, malcol)) +  guides(fill=guide_legend(title=expression(paste(italic("ndufs5")," Genotype"))))+
  labs(x = "Stages Behind Brood Maximum", y = "Count") +
  geom_histogram(bins = 6)
binned_chr13_relative
ggsave("output_files/CALL_malmito_chr13_embryos_relhist.pdf", binned_chr13_relative, 
       device = "pdf", width = 5, height = 2.92)


##########################
### Raw stage vs. ndufs5 figure (Sup.)
##########################

library(tidyverse)

mother_counts<-read.csv(file="input_files/CALL_mother_ndufs5_ndufa13_mito_genotypes.csv",
                        sep=",",head=TRUE, row.names = 1)

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

genotype_stage <- mother_embryo_stage_genotypes %>%
  mutate(mother = F) %>%
  rbind(mother_counts)

counts_by_stage <- genotype_stage %>%
  group_by(stage, mitotype, ndufs5) %>%
  summarize(Count = n()) %>%
  mutate(Stage = stage) %>%
  drop_na()
chr13_plot <- ggplot(counts_by_stage, aes(x = Stage, y = Count, fill = ndufs5)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.3),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(5.0, 7.5, 10.0, 12.0), labels = c("5.0","7.5","10.0","Maternal")) +
  scale_fill_manual(values = c(bircol, hetcol, malcol)) +
  guides(fill=guide_legend(title=expression(paste(italic("ndufs5")," Genotype"))))+
  facet_grid(mitotype ~ ., labeller = "label_both") + 
  geom_bar(stat = "identity")
chr13_plot
ggsave("output_files/CALL_chr13_embryos_raw.pdf", chr13_plot, 
       device = "pdf", width = 7, height = 6)


##########################
### Raw stage vs. ndufa13 figure (Sup.)
##########################

counts_by_stage <- genotype_stage %>%
  group_by(stage, mitotype, ndufa13) %>%
  summarize(Count = n()) %>%
  mutate(Stage = stage) %>%
  drop_na()
chr6_plot <- ggplot(counts_by_stage, aes(x = Stage, y = Count, fill = ndufa13)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.3),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(5.0, 7.5, 10.0, 12.0), labels = c("5.0","7.5","10.0","Maternal")) +
  scale_fill_manual(values = c(bircol, hetcol, malcol)) +
  guides(fill=guide_legend(title=expression(paste(italic("ndufa13")," Genotype"))))+
  facet_grid(mitotype ~ ., labeller = "label_both") + 
  geom_bar(stat = "identity")
chr6_plot
ggsave("output_files/CALL_chr6_embryos_raw.pdf", chr6_plot, 
       device = "pdf", width = 7, height = 6)


##########################
### Raw stage vs. Chr15 figure (Sup.)
##########################

counts_by_stage <- genotype_stage %>%
  group_by(stage, mitotype, ScyDAA6.5984.HRSCAF.6694.3367085) %>%
  summarize(Count = n()) %>%
  mutate(Stage = stage) %>%
  drop_na()
chr15_plot <- ggplot(counts_by_stage, aes(x = Stage, y = Count, fill = ScyDAA6.5984.HRSCAF.6694.3367085)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.3),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(5.0, 7.5, 10.0, 12.0), labels = c("5.0","7.5","10.0","Maternal")) +
  scale_fill_manual(values = c(bircol, hetcol, malcol)) +
  guides(fill=guide_legend(title="Chr 15 Genotype"))+
  facet_grid(mitotype ~ ., labeller = "label_both") + 
  geom_bar(stat = "identity")
chr15_plot

ggsave("output_files/CALL_chr15_embryos_raw.pdf", chr15_plot, device = "pdf", width = 7, height = 6)
