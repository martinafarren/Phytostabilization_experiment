#Biominerals Compost Pot Study
#qPCR statistics, Ending measurements

      # October 2019; Edits January 2023 & August 2023
      #MC Leewis mc.leewis@gmail.com

#load libraries
library(Rmisc)
library(ggplot2)
library(car)
library(rcompanion)
library("ggpubr")
library(ggplot2)
theme_set(theme_bw())

#Load Data
setwd("/Users/mleewis/OneDrive/Biominerals_Project/2019_CompostPotStudy/Endophytes/qPCR")
all.dat <- read.table("qPCR_tc.txt",row.names=1,header=T)

#subset data
compost <- subset(all.dat, Layer =="compost")
tailings <- subset(all.dat, Layer =="tailings")
initial <-  subset(all.dat, Layer =="initial")

######## COMPOST & TAILINGS comparisons #######
# These tests are all ONLY comparing the ending measurements. 
# All have un-even sample sizes, therefore can’t use parametric tests which assume even sample sizes

#KW for non-normal data, followed by wilcox post-hoc test
kruskal.test(Bact_gcn ~ trt, data = compost)
pairwise.wilcox.test(compost$Bact_gcn, compost$trt,
                     p.adjust.method = "fdr")

kruskal.test(Fung_gcn ~ trt, data = compost)
pairwise.wilcox.test(compost$Fung_gcn, compost$trt,
                     p.adjust.method = "fdr")

kruskal.test(FB_rati ~ trt, data = compost)
pairwise.wilcox.test(compost$FB_rati, compost$trt,
                     p.adjust.method = "fdr")

kruskal.test(Bact_gcn ~ trt, data = tailings)
pairwise.wilcox.test(tailings$Bact_gcn, tailings$trt,
                     p.adjust.method = "fdr")

kruskal.test(Fung_gcn ~ trt, data = tailings)
pairwise.wilcox.test(tailings$Fung_gcn, tailings$trt,
                     p.adjust.method = "fdr")

kruskal.test(FB_rati ~ trt, data = tailings)
pairwise.wilcox.test(tailings$FB_rati, tailings$trt,
                     p.adjust.method = "fdr")


##### qPCR figure ########
#1) Re-name the facet labels
layer.labs <- c("compost layer", "initial material", "tailings")
names(layer.labs) <- c("compost", "initial", "tailings")

#2) create each individual graph
b_p_comp  <- ggplot(compost, aes(x= trt, l.Bact_gcn)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = "Log Bacteria Gene Copy Number",
       subtitle = "(a) bacteria" #only compost gets the label, since it's first in line
  )+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(5,8.5)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))
b_p_comp

b_p_tailings <- ggplot(tailings, aes(x= trt, l.Bact_gcn)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = element_blank(), subtitle = " ")+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(5,8.5)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))
b_p_tailings

b_p_initial <- ggplot(initial, aes(x= trt, l.Bact_gcn)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = element_blank(), subtitle = " ")+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(5,8.5)+
  scale_x_discrete(labels = c("C(in)", "T(in)"))+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))
b_p_initial


f_p_comp  <- ggplot(compost, aes(x= trt, l.Fung_gcn)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = "Log Fungal Gene Copy Number",
       subtitle = "(b) fungi" #only compost gets the label, since it's first in line
  )+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(0,8.25)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))
f_p_comp

f_p_tailings <- ggplot(tailings, aes(x= trt, l.Fung_gcn)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = element_blank(), subtitle = " ")+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(0,8.25)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))

f_p_tailings

f_p_initial <- ggplot(initial, aes(x= trt, l.Fung_gcn)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = element_blank(), subtitle = " ")+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  scale_x_discrete(labels = c("C(in)", "T(in)"))+
  ylim(0,8.25)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))

f_p_initial

fb_p_comp  <- ggplot(compost, aes(x= trt, FB_rati)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = "Fungi : Bacteria", #only compost gets y-axis label, since it's first in line. 
       subtitle = "(c) fungi : bacteria" #only compost gets the label, since it's first in line
  )+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(0,2.0)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))

fb_p_comp

fb_p_tailings <- ggplot(tailings, aes(x= trt, FB_rati)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y=element_blank(), subtitle = " ")+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(0,0.65)+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))

fb_p_tailings


fb_p_initial <- ggplot(initial, aes(x= trt, FB_rati)) + 
  geom_boxplot(alpha = 0.1, aes(color = NULL)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .75)+ 
  theme_bw(base_size=13) +
  labs(x = "treatment", y = element_blank(), subtitle = " ")+ 
  theme(text = element_text(size = 10,face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10,  face = "plain"))+
  ylim(0,0.015)+
  scale_x_discrete(labels = c("C(in)", "T(in)"))+
  facet_grid(~Layer,
             labeller = labeller(Layer = layer.labs))

fb_p_initial


#3) Combine graphs into one figure 
  #layout: 
    #A# Bact_gcn: compost, tailings, initial
    #B# Fung_gcn: compost, tailings, initial
    #C# FB_rati : compost, tailings, initial

fig <- ggarrange(b_p_comp, b_p_tailings,  b_p_initial,
                 f_p_comp, f_p_tailings,  f_p_initial,
                 fb_p_comp, fb_p_tailings,  fb_p_initial,
                 ncol= 3, nrow = 3, 
                 widths  = c(0.29,0.45, 0.2))
fig

