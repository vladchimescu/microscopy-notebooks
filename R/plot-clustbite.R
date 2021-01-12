## Plot heatmap for BiTE
library(pheatmap)
library(dplyr)
library(grid)

datadir = "~/Documents/github/microscopy-notebooks/data/"
figdir = "~/Documents/github/microscopy-notebooks/figures/"

save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

drug_prof = read.csv(paste0(datadir, 'all_clustbite_profiles.csv'),
                     row.names = 1)
drug_prof = drug_prof[rowSums(is.na(drug_prof)) < 500,]

# remove Hoechst
drug_prof = drug_prof[,!grepl("Hoechst", colnames(drug_prof))]
# remove filled_area
drug_prof = drug_prof[,!grepl("filled_area", colnames(drug_prof))]
# remove extent
drug_prof = drug_prof[,!grepl("extent", colnames(drug_prof))]
# remove solidity
drug_prof = drug_prof[,!grepl("solidity", colnames(drug_prof))]

drug_prof = drug_prof[!rownames(drug_prof)%in% c("DMSO"),]

colannot = data.frame(names = colnames(drug_prof))

colannot = tidyr::separate(colannot, names, 
                           into=c('tmp', 'Channel', 
                                  'feature', 'PatientID'),
                           sep = '\\.', 
                           remove=F)
colannot = select(colannot, names, Channel, feature, PatientID)

bitedir = "~/Documents/embl/gitlab/microscopy/data/Tobias/"
samples = readxl::read_xlsx(paste0(bitedir, "SampleList.xlsx"))

colannot = inner_join(colannot,
                      select(samples, PatientID, Entity))

rownames(colannot) = colannot$names
colannot = select(colannot, -c(names, PatientID))

chan = c( "#d8345f", "#63b7af", "#ffc93c")
names(chan) = c("APC", "Calcein", "PE")
dg = c("#010059", "#fdbfb3", "#fcf594", "#aee7e8")
names(dg) = c("CLL", "DLBCL", "FL", "MCL")
ft = c("#e7dfd5", "#84a9ac")
names(ft) = c("mean_intensity", "area")

ann_colors <- list(Channel=chan,
                   Entity=dg,
                   feature=ft)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#003836", "white", "#a87b00"))(paletteLength)
maxval = 5
myBreaks <- c(seq(-maxval, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, maxval,
                  length.out=floor(paletteLength/2)))
row_clust = pheatmap(cor(t(drug_prof), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(drug_prof, use="pairwise.complete.obs"), silent = T)

# plot the big heatmap
ph = pheatmap(drug_prof,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'clust-BiTE-drugprofiles-all.pdf'),
                  width = 9, height = 12)


drug_prof = drug_prof[rowSums(is.na(drug_prof)) < 75,]
row_clust = pheatmap(cor(t(drug_prof), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(drug_prof, use="pairwise.complete.obs"), silent = T)

# plot the heatmap for complete observations only
ph = pheatmap(drug_prof,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'clust-BiTE-drugprofiles-notmissing.pdf'),
                  width = 9, height = 7)

drug_prof = read.csv(paste0(datadir, 'clustbite_comb_profiles.csv'),
                     row.names = 1)
# remove Hoechst
drug_prof = drug_prof[,!grepl("Hoechst", colnames(drug_prof))]
# remove filled_area
drug_prof = drug_prof[,!grepl("filled_area", colnames(drug_prof))]
# remove extent
drug_prof = drug_prof[,!grepl("extent", colnames(drug_prof))]
# remove solidity
drug_prof = drug_prof[,!grepl("solidity", colnames(drug_prof))]

drug_prof = drug_prof[!rownames(drug_prof)%in% c("DMSO"),]
row_clust = pheatmap(cor(t(drug_prof), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(drug_prof, use="pairwise.complete.obs"), silent = T)

# plot the heatmap for complete observations only
ph = pheatmap(drug_prof,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'clust-BiTE-combprofiles-all.pdf'),
                  width = 9, height = 8)

drug_prof = drug_prof[rowSums(is.na(drug_prof)) < 75,]
row_clust = pheatmap(cor(t(drug_prof), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(drug_prof, use="pairwise.complete.obs"), silent = T)

# plot the heatmap for complete observations only
ph = pheatmap(drug_prof,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'clust-BiTE-combprofiles-notmissing.pdf'),
                  width = 9, height = 5)

df = read.csv(paste0(datadir, 'BiTE/clustbite_all.csv'))

topwells = select(df, comb, PatientID, 
       ch.Calcein.area, 
       ch.Calcein.mean_intensity,
       ch.PE.area,
       ch.PE.mean_intensity,
       ch.APC.area,
       ch.APC.mean_intensity) %>%
top_n(n=30, ch.Calcein.area) %>%
  arrange(desc(ch.Calcein.area))

write.csv(topwells, 
          file = paste0(datadir, 'bite-topwells.csv'),
          quote = F, row.names = F)

drug_prof$comb = rownames(drug_prof)
df = reshape2::melt(drug_prof, id.vars="comb") %>%
  tidyr::separate(variable, into = c("tmp", "channel", "feature", "PatientID"), sep = "\\.") %>%
  select(-tmp)

library(ggplot2)
p = ggplot(filter(df, comb %in% c("CD19", "CD20", 
                              "BCMA", "CD19_Lenalidomide",
                              "CD20_Lenalidomide")),
       aes(x = comb, y = value, fill = comb))+
  geom_boxplot(alpha=0.5) +
  geom_dotplot(stackdir = "center", binaxis = 'y') + 
  labs(x = "", y = "Normalized values") + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  facet_wrap(~ channel + feature, ncol = 2) +
  theme_bw() +
  coord_fixed() + 
  theme(axis.text.x = element_text(angle=45,hjust=1),
        aspect.ratio = 1,
        legend.position = 'none' )
ggsave(p, filename = paste0(figdir, "bite-boxplots.pdf"),
       height = 9, width = 7)
