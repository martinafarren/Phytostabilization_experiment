#Biominerals Compost Pot Study
# Microbiome data analysis

## Permanova, 
## Pairwise Permanova, 
## Differential expression analysis (DESeq2), 
## Plot of Relative abundance of bacterial (a) and fungal (b) phyla 

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
output <- "c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/finalni_vystup"
setwd(output)

# IMPORTANT: DATA ANALYSES ARE SHOWN ON AN EXAMPLE OF TAILINGS SAMPLES, specifically bacterial communities
# Same analyses were run also on samples of compost and roots

#### LOAD Data AFTER rarefy ####
#  Tailings: 16S
ps_tailings_16S <- readRDS("c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input/ps_tailings_16S.rds")


#### LOAD Data BEFORE rarefy ####
#  Tailings: 16S
ps_tailings_16S_beforerarefy <- readRDS("c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input/ps_tailings_16S_BEFORErarefy.rds")

#### Load information about samples ####
ps_USGS_16S_Biominerals <- readRDS("c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input/phyloseq_USGS_16S_Biominerals.rds")
env.data_16S <- read.csv("c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input/env.data_16s_Silva.csv", header = T, row.names=1)
env.data_16S_biominerals <- env.data_16S[env.data_16S$Project=="Biominerals_project",]
sample_data(ps_USGS_16S_Biominerals)$treatment <- factor(sample_data(ps_USGS_16S_Biominerals)$treatment_II, 
                                                         levels = c("A_initial material", "B_T", "C_TP", "D_TPE", "E_TC", "F_TPC", "G_TPEC"), 
                                                         labels = c("initial material", "T", "TP", "TPE", "TC", "TPC", "TPEC"))
env_tailings_initials_16S <- env.data_16S_biominerals[env.data_16S_biominerals$material=="tailing",]
env_tailings_16S <- env_tailings_initials_16S[env_tailings_initials_16S$treatment !="initial material",]
rownames(env_tailings_16S)==sample_names(ps_tailings_16S)


#### DATA ANALYSES ####
#### PERMANOVA ####
# 16S
Hellinger_tailings_16S = decostand(otu_table(ps_tailings_16S), "hell")        		#Hellinger transform
s.dis_Hellinger_tailings_16S = vegdist(Hellinger_tailings_16S , "bray") 	   		#Bray-Curtis distance
adonis(formula =s.dis_Hellinger_tailings_16S~env_tailings_16S$treatment, permutations = 999)                            
adonis(formula =s.dis_Hellinger_tailings_16S~env_tailings_16S$coating, permutations = 999)                  
adonis(formula =s.dis_Hellinger_tailings_16S~env_tailings_16S$compost, permutations = 999)            
adonis(formula =s.dis_Hellinger_tailings_16S~env_tailings_16S$plant, permutations = 999)              

length(env_tailings_16S$treatment)==length(rownames(Hellinger_tailings_16S))

adonis(formula =s.dis_Hellinger_tailings_16S~compost*coating*plant, data = env_tailings_16S, permutations = 999)                             # 0.001 ***

bd_l <- betadisper(s.dis_Hellinger_tailings_16S,env_tailings_16S$treatment)
permutest(bd_l)

#### post hoc adonis ####
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~ 
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))], permutations = 9999, method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
} 

# Now run the multiple comparisons - be aware of using full OTU matrix, not distance matrix in this step
pairwise.adonis(Hellinger_tailings_16S,env_tailings_16S$treatment, sim.method="bray",p.adjust.m = "fdr")



##### DeSeq ######
## Run on data before rarefy
# 16S
taxa_names(ps_tailings_16S_beforerarefy)
taxa_names(ps_tailings_16S_beforerarefy)<- paste(taxa_names(ps_tailings_16S_beforerarefy), "_",tax_table(ps_tailings_16S_beforerarefy)[,6])
taxa.rename <- as.data.frame(tax_table(ps_tailings_16S_beforerarefy))


## Agglomerate ASVs to Genus level ##
ps <- tax_glom(ps_tailings_16S_beforerarefy, "Genus")
sample_data(ps)
length(get_taxa_unique(ps, taxonomic.rank = "Genus"))

# Model
dds <- phyloseq_to_deseq2(ps, ~ treatment)
dds <- dds[rowSums(counts(dds)) > 1, ]

# estimation of size factors, estimation of dispersion, Negative Binomial GLM fitting and Wald statistics
dds <- DESeq(dds)

# results
res = results(dds, cooksCutoff = FALSE)
res

# lfcShrink applied to shrink logarithmic fold change values
# example: T vs TPC conrast
sig_shrink_T_vs_TPC <- lfcShrink(dds, contrast = c("treatment", "T", "TPC"), type="normal")
sig_shrink_T_vs_TPC


# only significant results
alpha = 0.01
sig_shrink_T_vs_TPC = sig_shrink_T_vs_TPC[order(sig_shrink_T_vs_TPC$padj, na.last=NA), ]
sigtab = sig_shrink_T_vs_TPC[(sig_shrink_T_vs_TPC$padj < alpha), ]
# assign the taxonomy
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))

# create a table with results
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
write.csv(sigtabgen, "Deseq_T_vs_TPC.csv")




#### BARPLOT PHYLA ####
#  1) TAXGLOM
ps_16Sall_Phylum <- tax_glom(ps_USGS_16S_Biominerals, "Phylum")


#  2) COLORS
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
PhylumList = unique(tax_table(ps_16Sall_Phylum)[,"Phylum"])
PhylumPalette = getPalette(length(PhylumList))
names(PhylumPalette) = PhylumList

#  3)  tailingS - MERGING/TRANSFORMATION
merge_16S_tailing_Phylum <- merge_samples(ps_16S_tailing, sample_data(ps_16S_tailing)$unique, fun = mean)
sample_data(merge_16S_tailing_Phylum)
transformed_counts_tailing <- transform_sample_counts(merge_16S_tailing_Phylum, function(x){x / sum(x)})
# "Other" for < 0.1%
tax_table(transformed_counts_tailing)[taxa_sums(transformed_counts_tailing)/sum(taxa_sums(transformed_counts_tailing))<=0.001,2]<- "Other"
env_tailing <- data.frame(treatment = c("B_T","E_TC","C_TP","F_TPC","D_TPE","G_TPEC"),
                          material = rep("tailings", times=6))
rownames(env_tailing)  <- sample_names(transformed_counts_tailing)
sample_data(transformed_counts_tailing) <- env_tailing
sample_data(transformed_counts_tailing)$treatment <- factor(sample_data(transformed_counts_tailing)$treatment, 
                                                            levels = c("B_T", "C_TP", "D_TPE", "E_TC", "F_TPC", "G_TPEC"), 
                                                            labels = c("T", "TP", "TPE", "TC", "TPC", "TPEC"))

#  4)  tailing - BARPLOT
Plot_16S_tailing <- plot_bar(transformed_counts_tailing, "treatment",fill = "Phylum", facet_grid = ~material)+
  theme_bw(base_size = 13)+
  labs(subtitle = " ")+
  #scale_fill_brewer(palette = "Set2")+
  scale_fill_manual(values = PhylumPalette) +
  theme(text = element_text(size = 10,face = "plain"), 
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
Plot_16S_tailing

