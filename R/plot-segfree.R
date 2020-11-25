## Plots for organoid screen
library(dplyr)
library(ggplot2)
library(pheatmap)
datadir = '~/Documents/github/microscopy-notebooks/data/segfree/'

df = read.csv(paste0(datadir, 'borten-brc-profiles.csv'), row.names = 1)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#003836", "white", "#a87b00"))(paletteLength)
maxval = 5
myBreaks <- c(seq(-maxval, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, maxval,
                  length.out=floor(paletteLength/2)))
row_clust = pheatmap(cor(t(df), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(df, use="pairwise.complete.obs"), silent = T)

ph = pheatmap(df,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              treeheight_row = 0, 
              treeheight_col = 0)

#df = df[,grepl('superblock',colnames(df))]
