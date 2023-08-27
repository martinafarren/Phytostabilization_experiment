#Biominerals Compost Pot Study
# Microbial diversity Plots and Statistics

# October 2019; Edits August 2023
# Martina Farren (Kracmarova), kracma.mk@gmail.com

rm(list=ls())
#install.packages("tidyverse")
#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("dada2","phyloseq","phangorn"))
library(phyloseq); packageVersion("phyloseq")
library(vegan)
library(devEMF)
library(ggplot2)
require(gridExtra)
library(DESeq2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
set.seed(1)


#### LOADING DATA####
output <- "c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/finalni_vystup/level2"
setwd(output)

ps_16S_Biominerals_beforerarefy <- readRDS("c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input/ps_USGS_16S_Biominerals_predrarefy.rds")
env.data_16S <- read.csv("c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input/env.data_16s_Silva.csv", header = T, row.names=1)
sample_data(ps_16S_Biominerals_beforerarefy)$treatment <- factor(sample_data(ps_16S_Biominerals_beforerarefy)$treatment_II, 
                                                         levels = c("A_initial material", "B_T", "C_TP", "D_TPE", "E_TC", "F_TPC", "G_TPEC"), 
                                                         labels = c("initial material", "T", "TP", "TPE", "TC", "TPC", "TPEC"))
sample_data(ps_16S_Biominerals_beforerarefy)$material <- factor(sample_data(ps_16S_Biominerals_beforerarefy)$material, 
                                                        levels = c("compost", "root", "tailing"), 
                                                        labels = c("compost layer", "plant roots", "tailings"))


ps_ITS_Biominerals_beforerarefy <- readRDS("C:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input_ITS/ps_USGS_ITS_Biominerals_predrarefy.rds")
sample_data(ps_ITS_Biominerals_beforerarefy)
env.data_ITS <- read.csv("C:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input_ITS/env_USGS_ITS.csv", header = T, row.names=1)
sample_data(ps_ITS_Biominerals_beforerarefy)$treatment <- factor(sample_data(ps_ITS_Biominerals_beforerarefy)$treatment_II, 
                                                         levels = c("A_initial material", "B_T", "C_TP", "D_TPE", "E_TC", "F_TPC", "G_TPEC"), 
                                                         labels = c("initial material", "T", "TP", "TPE", "TC", "TPC", "TPEC"))

sample_data(ps_ITS_Biominerals_beforerarefy)$material <- factor(sample_data(ps_ITS_Biominerals_beforerarefy)$material, 
                                                        levels = c("compost", "root", "tailing"), 
                                                        labels = c("compost layer", "plant roots", "tailings"))



#### SUBSETS OF SAMPLES ####
#### TAILINGS, COMPOST, ROOTS, INITIAL MATERIALs ####
ps_USGS_16S_Biominerals = ps_16S_Biominerals_beforerarefy
ps_USGS_ITS_Biominerals = ps_ITS_Biominerals_beforerarefy

# 16S - tailings, compost, root
ps_16S_without_initials <- subset_samples(ps_USGS_16S_Biominerals, treatment !="initial material")
ps_16S_without_initials <- prune_taxa(taxa_sums(ps_16S_without_initials)>0, ps_16S_without_initials)
ps_16S_without_initials
sample_data(ps_16S_without_initials)

ps_16S_tailing <- subset_samples(ps_16S_without_initials, material =="tailings")
ps_16S_compost <- subset_samples(ps_16S_without_initials, material =="compost layer")
ps_16S_root <- subset_samples(ps_16S_without_initials, material =="plant roots")

###
# ITS - tailings, compost, root
ps_ITS_without_initials <- subset_samples(ps_USGS_ITS_Biominerals, treatment !="initial material")
ps_ITS_without_initials <- prune_taxa(taxa_sums(ps_ITS_without_initials)>0, ps_ITS_without_initials)
ps_ITS_without_initials
sample_data(ps_ITS_without_initials)

ps_ITS_tailing <- subset_samples(ps_ITS_without_initials, material =="tailings")
ps_ITS_compost <- subset_samples(ps_ITS_without_initials, material =="compost layer")
ps_ITS_root <- subset_samples(ps_ITS_without_initials, material =="plant roots")


#### SHANNON BACTERIAL COMMUNITIES ####

# Shannon - 16S - tailing
(p_16S_tailing = plot_richness(ps_16S_tailing, x = c("treatment"),measures = c("Shannon")))                
pp_16S_tailing = p_16S_tailing + geom_boxplot(data = p_16S_tailing$data, aes(x = treatment, y = value, color = NULL), alpha = 0.1)  
pp_16S_tailing <- pp_16S_tailing + 
  facet_grid(~material)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(2, 5)+
  #labs(title = "Shannon diversity index", subtitle = "(a) bacteria")+
  labs(subtitle = " ")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

pp_16S_tailing

# Shannon - 16S - compost
(p_16S_compost = plot_richness(ps_16S_compost, x = c("treatment"),measures = c("Shannon")))                
pp_16S_compost = p_16S_compost + geom_boxplot(data = p_16S_compost$data, aes(x = treatment, y = value, color = NULL), alpha = 0.1)  
pp_16S_compost <- pp_16S_compost + 
  facet_grid(~material)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(2, 5)+
  #labs(title = "Shannon diversity index", subtitle = "(a) bacteria")+
  labs(subtitle = "(a) bacteria")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"), # angle = 90
        legend.text = element_text(size = 10,  face = "plain"))
pp_16S_compost

# Shannon - 16S - root
(p_16S_root = plot_richness(ps_16S_root, x = c("treatment"),measures = c("Shannon")))                
pp_16S_root = p_16S_root + geom_boxplot(data = p_16S_root$data, aes(x = treatment, y = value, color = NULL), alpha = 0.1)  
pp_16S_root <- pp_16S_root + 
  facet_grid(~material)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(2, 5)+
  #labs(title = "Shannon diversity index", subtitle = "(a) bacteria")+
  labs(subtitle = " ")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
pp_16S_root


#### SHANNON FUNGAL COMMINITIES ####
# Shannon - ITS - tailing
(p_ITS_tailing = plot_richness(ps_ITS_tailing, x = c("treatment"),measures = c("Shannon")))                
pp_ITS_tailing = p_ITS_tailing + geom_boxplot(data = p_ITS_tailing$data, aes(x = treatment, y = value, color = NULL), alpha = 0.1)  
pp_ITS_tailing <- pp_ITS_tailing + 
  facet_grid(~material)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(0.7, 4)+
  #labs(title = "Shannon diversity index", subtitle = "(a) bacteria")+
  labs(subtitle = " ")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

pp_ITS_tailing

# Shannon - ITS - compost
(p_ITS_compost = plot_richness(ps_ITS_compost, x = c("treatment"),measures = c("Shannon")))                
pp_ITS_compost = p_ITS_compost + geom_boxplot(data = p_ITS_compost$data, aes(x = treatment, y = value, color = NULL), alpha = 0.1)  
pp_ITS_compost <- pp_ITS_compost + 
  facet_grid(~material)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(0.7, 4)+
  #labs(title = "Shannon diversity index", subtitle = "(a) bacteria")+
  labs(subtitle = "(b) fungi")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))
pp_ITS_compost

# Shannon - ITS - root
(p_ITS_root = plot_richness(ps_ITS_root, x = c("treatment"),measures = c("Shannon")))                
pp_ITS_root = p_ITS_root + geom_boxplot(data = p_ITS_root$data, aes(x = treatment, y = value, color = NULL), alpha = 0.1)  
pp_ITS_root <- pp_ITS_root + 
  facet_grid(~material)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(0.7, 4)+
  #labs(title = "Shannon diversity index", subtitle = "(a) bacteria")+
  labs(subtitle = " ")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
pp_ITS_root

#### SHANNON INITIAL MATERIALS: bacterial and fungal ####
# 16S
ps_16S_initials <- subset_samples(ps_USGS_16S_Biominerals, treatment =="initial material")
sample_data(ps_16S_initials)
env.data_16S_initials <- env.data_16S[env.data_16S$treatment=="initial material",]
env.data_16S_initials
env.data_16S_initials <- env.data_16S_initials %>% 
  unite(material_initials, c("material", "coating"))
sample_data(ps_16S_initials)$material_initials <- as.factor(env.data_16S_initials$material_initials)
sample_data(ps_16S_initials)$material_initials <- factor(sample_data(ps_16S_initials)$material_initials, 
                                                         levels = c("compost_not_defined", "root_added_endophytes", "root_NO_added_endophytes",
                                                                    "tailing_not_defined"), 
                                                         labels = c("C(in)", "S(+e)", "S(-e)","T(in)"))

# Shannon - 16S
(p_16S_initials = plot_richness(ps_16S_initials, x = c("material_initials"),measures = c("Shannon")))                
pp_16S_initials = p_16S_initials + geom_boxplot(data = p_16S_initials$data, aes(x = material_initials, y = value, color = NULL), alpha = 0.1)  
plot_16S_initials <- pp_16S_initials + 
  facet_grid(~treatment)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(2, 5)+
  labs(subtitle = " ")+
  xlab("initial material")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
plot_16S_initials

# ITS
ps_ITS_initials <- subset_samples(ps_USGS_ITS_Biominerals, treatment =="initial material")
sample_data(ps_ITS_initials)
env.data_ITS_initials <- env.data_ITS[env.data_ITS$treatment=="initial material",]
env.data_ITS_initials
env.data_ITS_initials <- env.data_ITS_initials %>% 
  unite(material_initials, c("material", "coating"))
sample_data(ps_ITS_initials)$material_initials <- as.factor(env.data_ITS_initials$material_initials)
sample_data(ps_ITS_initials)$material_initials <- factor(sample_data(ps_ITS_initials)$material_initials, 
                                                         levels = c("compost_not_defined", "root_added_endophytes", "root_NO_added_endophytes",
                                                                    "tailing_not_defined"), 
                                                         labels = c("C(in)", "S(+e)", "S(-e)","T(in)"))

# Shannon - ITS
(p_ITS_initials = plot_richness(ps_ITS_initials, x = c("material_initials"),measures = c("Shannon")))                
pp_ITS_initials = p_ITS_initials + geom_boxplot(data = p_ITS_initials$data, aes(x = material_initials, y = value, color = NULL), alpha = 0.1)  
plot_ITS_initials <- pp_ITS_initials + 
  facet_grid(~treatment)+
  theme_bw(base_size = 13)+
  #geom_point(aes(color=coating))+
  ylim(0.7, 4)+
  labs(subtitle = " ")+
  xlab("initial material")+
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
plot_ITS_initials


#### FINAL PLOTs ####
grid.arrange(pp_16S_compost, pp_16S_root,  pp_16S_tailing, plot_16S_initials,
             pp_ITS_compost, pp_ITS_root, pp_ITS_tailing,  plot_ITS_initials, 
             ncol= 4, widths  = c(0.29,0.17,0.45, 0.30))



#### Pairwise Wilcox test: TAILINGS with initial material ####
library(multcomp)
library(FSA)
# 16S
sample_data(ps_USGS_16S_Biominerals)
ps_tailings_with_init_16S <- subset_samples(ps_USGS_16S_Biominerals, material == "tailings")

richness.values_tailings_16S <- estimate_richness(ps_tailings_with_init_16S, measures = c("Shannon"))
pairwise.wilcox.test(richness.values_tailings_16S$Shannon, sample_data(ps_tailings_with_init_16S)$treatment,
                     p.adjust.method = "fdr")


# ITS
sample_data(ps_USGS_ITS_Biominerals)
ps_tailings_with_init_ITS <- subset_samples(ps_USGS_ITS_Biominerals, material == "tailings")

richness.values_tailings_ITS <- estimate_richness(ps_tailings_with_init_ITS, measures = c("Shannon"))
pairwise.wilcox.test(richness.values_tailings_ITS$Shannon, sample_data(ps_tailings_with_init_ITS)$treatment,
                     p.adjust.method = "fdr")


#### Pairwise Wilcox test: COMPOST with initial material ####

# 16S
ps_compost_with_init_16S <- subset_samples(ps_USGS_16S_Biominerals, material == "compost layer")
sample_data(ps_compost_with_init_16S)

richness.values_compost_ITS <- estimate_richness(ps_compost_with_init_16S, measures = c("Shannon"))
pairwise.wilcox.test(richness.values_compost_ITS$Shannon, sample_data(ps_compost_with_init_16S)$treatment,
                     p.adjust.method = "fdr")

# ITS
ps_compost_with_init_ITS <- subset_samples(ps_USGS_ITS_Biominerals, material == "compost layer")
sample_data(ps_compost_with_init_ITS)

richness.values_compost_ITS <- estimate_richness(ps_compost_with_init_ITS, measures = c("Shannon"))
pairwise.wilcox.test(richness.values_compost_ITS$Shannon, sample_data(ps_compost_with_init_ITS)$treatment,
                     p.adjust.method = "fdr")

#### Pairwise Wilcox test: ROOTS with initial material ####

# 16S
ps_root_with_init_16S <- subset_samples(ps_USGS_16S_Biominerals, material == "plant roots")
sample_data(ps_root_with_init_16S)

richness.values_root_ITS <- estimate_richness(ps_root_with_init_16S, measures = c("Shannon"))
pairwise.wilcox.test(richness.values_root_ITS$Shannon, sample_data(ps_root_with_init_16S)$treatment,
                     p.adjust.method = "fdr")

# ITS
ps_root_with_init_ITS <- subset_samples(ps_USGS_ITS_Biominerals, material == "plant roots")
sample_data(ps_root_with_init_ITS)

richness.values_root_ITS <- estimate_richness(ps_root_with_init_ITS, measures = c("Shannon"))
pairwise.wilcox.test(richness.values_root_ITS$Shannon, sample_data(ps_root_with_init_ITS)$treatment,
                     p.adjust.method = "fdr")






