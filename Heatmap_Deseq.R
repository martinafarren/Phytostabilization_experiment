#### Heatmap based on DeSeq results ####


rm(list=ls())

#install.packages("tidyverse")
#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("dada2","phyloseq","phangorn"))
library(vegan)
library(devEMF)
library(ggplot2)
library(extrafont)
library("gplots")
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
set.seed(1)


path <- "c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/input_heatmap"
output <- "c:/Users/kracmarm/R_analyzy/Kracma/USGS_Biominerals/finalni_vystup"
setwd(output)

## File import POSITIVE and NEGATIVE VALUES (LOGFOL2)##
data_allVal <- read.csv(file.path(path, "DeSeq_16S_allValues_lcfShrink.csv"), header = T, row.names = 1)
data_allVal
dim(data_allVal)
data_allVal[duplicated(data_allVal$X), ]
any(duplicated(data_allVal$X))
rownames(data_allVal)

data_allVal.matrx <- as.matrix(data_allVal)
data_allVal.matrx[is.na(data_allVal.matrx)] <- 0
class(data_allVal.matrx)


# changes the distance measure and clustering method
row_distance = dist(data_allVal.matrx, method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(data_allVal.matrx), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")

# names in italic
make_italic_names <- function(mat, rc_fun, rc_names) {
  italic_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        italic_names[i] <<-
        bquote(italic(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  italic_names
}

# heatmap
pheatmap(data_allVal.matrx,
         labels_row = make_italic_names(data_allVal.matrx, rownames, row.names(data_allVal)),
         cexRow = 0.5,                # size of text in rows
         cexCol = 1,                # size of text in columns
         angle_col = 0,
         cluster_rows = F,                  # NO dendogram on rows
         clustering_distance_cols = "euclidean",   #  dendogram on columns: manhattan, canberra, 
         clustering_method = "ward.D2",    #  ward.D2, 
         #cellnote = data_allVal.matrx,  # same data_allVal set for cell labels
         #main = "   ",             # heat map title
         fontsize_row = 10,         # font size for row labels
         fontsize_col = 11,         # font size for column labels
         #cellheight = 8,
         #cellwidth = NA,
         #density.info = "none",  # turns off density plot inside color legend
         border_color=NA,         # turns off trace lines inside the heat map
         na_col = "white",
         #margins = c(5,7),     # widens margins around plot
         #gaps_row = 20,         # place gap under x-th row
         #color = redblue(100),
)




