## Plot the heatmap of all features
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

drug_prof = read.csv(paste0(datadir, 'all_profiles_coculture.csv'),
                     row.names = 1)

colannot = data.frame(names = colnames(drug_prof))
colannot = mutate(colannot, Channel = ifelse(grepl('Calcein', names), 'Calcein',
                                             ifelse(grepl('Hoechst', names), 'Hoechst', 'Lysosomal')))
colannot = mutate(colannot, Culture = ifelse(grepl("_C\\.", names), "Coculture", "Monoculture"))
colannot = mutate(colannot, plate = gsub("(.+_[CM]\\.)(18.+)", '\\2', names))

# add diagnosis based on the screened plate
patannot = readxl::read_xlsx("~/Documents/embl/gitlab/microscopy/data/Metafiles/Screened_Patients.xlsx")
patannot = filter(patannot, Diagnosis != "CLL") %>%
  mutate(Screen = gsub("181106", "181109", Screen),
         Plate_Nr = gsub("_", "", Plate_Nr)) %>%
  mutate(plate = paste(Screen, Plate_Nr, sep = "_")) %>%
  select(plate, Diagnosis)
colannot = left_join(colannot, patannot)

rownames(colannot) = colannot$names
colannot = select(colannot, -c(names, plate))

chan = c("#63b7af", "#035aa6", "#d8345f")
names(chan) = c("Calcein", "Hoechst", "Lysosomal")
dg = c("#010059", "#fdbfb3", "#fcf594", "#aee7e8")
names(dg) = c("AML", "HCL", "MCL", "T-PLL")
# for culture
colscult = c("#407088", "#b3b3b3")
names(colscult) = c("Coculture", "Monoculture")
ann_colors <- list(Channel=chan,
                   Diagnosis=dg,
                   Culture=colscult)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#003836", "white", "#a87b00"))(paletteLength)
maxval = 6
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

save_pheatmap_pdf(ph, filename = paste0(figdir, 'drugprofiles-all.pdf'),
                  width = 9, height = 20)


# now plot the heatmap for a selection of drugs
drug_sel = c('Tofacitinib', 'Midostaurin',
            'Ganetespib', 'Lenalidomide',
            'Pyridone 6', 'UMI-77',
            'Bafilomycin A1',
            'Quizartinib', 'Hydroxychloroquine',
            'Fludarabine', 'Vorinostat',
            'Thioguanine', 'Nutlin 3a',
            'Palbociclib', 'Carfilzomib',
            'JQ1', 'Cytarabine',
            'BAY61-3606', 'Venetoclax',
            'Ixazomib')

drugprof_sel = drug_prof[gsub("(.+)_(.+)", "\\1", rownames(drug_prof)) %in% drug_sel,]
drugprof_sel = drugprof_sel[rowSums(is.na(drugprof_sel)) < 1000,]
row_clust = pheatmap(cor(t(drugprof_sel), use="pairwise.complete.obs"), silent = T)
#col_clust = pheatmap(cor(drugprof_sel, use="pairwise.complete.obs"), silent = T)

# plot the big heatmap
ph = pheatmap(drugprof_sel,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'drugprofiles-selected.pdf'),
                  width = 9, height = 10)

# even smaller subset for the presentation
drug_sel = c('Tofacitinib', 'Midostaurin',
             'Ganetespib', 'Lenalidomide',
             'Pyridone 6', 'UMI-77',
             'Bafilomycin A1',
             'Quizartinib', 'Hydroxychloroquine',
             'Vorinostat',
             'Thioguanine', 'Nutlin 3a',
             'Palbociclib', 'Carfilzomib',
             'JQ1', 'Venetoclax')

drugprof_sel = drug_prof[gsub("(.+)_(.+)", "\\1", rownames(drug_prof)) %in% drug_sel,]
drugprof_sel = drugprof_sel[rowSums(is.na(drugprof_sel)) < 1000,]
row_clust = pheatmap(cor(t(drugprof_sel), use="pairwise.complete.obs"), silent = T)
#col_clust = pheatmap(cor(drugprof_sel, use="pairwise.complete.obs"), silent = T)

# plot the big heatmap
ph = pheatmap(drugprof_sel,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'drugprofiles-small-subset.pdf'),
                  width = 9, height = 8)

# Coculture vs Monoculture morphological feature difference by sample
cultdiff = read.csv(paste0(datadir, 'diff_DMSO.csv'))
cultdiff = dplyr::rename(cultdiff, feature=X)
cultdiff = mutate(cultdiff, diff_medians = ifelse(padj < 0.01, diff_medians, 0))

df_wide = tidyr::spread(select(cultdiff, feature, plate, diff_medians),
                        key='feature', value='diff_medians')
rownames(df_wide) = df_wide$plate
df_wide = df_wide[,-1]
colnames(df_wide) = gsub("ch-", "", colnames(df_wide))

df_wide = t(df_wide)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#78A1A7","white", "#C9788A"))(paletteLength)
maxval = max(abs(df_wide), na.rm = T)
myBreaks <- c(seq(-maxval, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, maxval,
                  length.out=floor(paletteLength/2)))

row_clust = pheatmap(cor(t(df_wide), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(df_wide, use="pairwise.complete.obs"), silent = T)


colannot = data.frame(plate = colnames(df_wide))
colannot = left_join(colannot, patannot)
rownames(colannot) = colannot$plate
colannot = select(colannot, -c(plate))
ann_colors <- list(Diagnosis=dg)

df_wide[df_wide == 0] = NA

ph = pheatmap(df_wide,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0,
              na_col = '#888888')
save_pheatmap_pdf(ph, filename = paste0(figdir, 'DMSO-culture-difference.pdf'),
                  width = 8, height = 2)


# autophagy drugs
drug_sel = c('Midostaurin', 'Ganetespib', 
             'Lenalidomide', 'UMI-77',
             'Bafilomycin A1', 'Quizartinib', 
             'Hydroxychloroquine','Vorinostat',
             'Carfilzomib','JQ1', 
             'Venetoclax', 'Ixazomib',
             'MK2206', 'PRT062607',
             'IRAK4 Inhibitor, Compound 26')

drugprof_sel = drug_prof[gsub("(.+)_(.+)", "\\1", rownames(drug_prof)) %in% drug_sel,]
drugprof_sel = drugprof_sel[rowSums(is.na(drugprof_sel)) < 1000,]
row_clust = pheatmap(cor(t(drugprof_sel), use="pairwise.complete.obs"), silent = T)
#col_clust = pheatmap(cor(drugprof_sel, use="pairwise.complete.obs"), silent = T)

# plot the big heatmap
ph = pheatmap(drugprof_sel,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'drugprofiles-autophagy.pdf'),
                  width = 9, height = 8)

cultdiff = read.csv(paste0(datadir, 'diff_DMSO.csv'))
cultdiff = dplyr::rename(cultdiff, feature=X)
cultdiff = mutate(cultdiff, diff_medians = ifelse(padj < 0.01, diff_medians, 0))

df_wide = tidyr::spread(select(cultdiff, feature, plate, diff_medians),
                        key='feature', value='diff_medians')
rownames(df_wide) = df_wide$plate
df_wide = df_wide[,-1]
colnames(df_wide) = gsub("ch-", "", colnames(df_wide))

df_wide = t(df_wide)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#78A1A7","white", "#C9788A"))(paletteLength)
maxval = max(abs(df_wide), na.rm = T)
myBreaks <- c(seq(-maxval, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, maxval,
                  length.out=floor(paletteLength/2)))

row_clust = pheatmap(cor(t(df_wide), use="pairwise.complete.obs"), silent = T)
col_clust = pheatmap(cor(df_wide, use="pairwise.complete.obs"), silent = T)


colannot = data.frame(plate = colnames(df_wide))
colannot = left_join(colannot, patannot)
rownames(colannot) = colannot$plate
colannot = select(colannot, -c(plate))
ann_colors <- list(Diagnosis=dg)

df_wide[df_wide == 0] = NA

ph = pheatmap(df_wide,
              cluster_rows = row_clust$tree_row, 
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0,
              na_col = '#888888')
save_pheatmap_pdf(ph, filename = paste0(figdir, 'DMSO-culture-difference.pdf'),
                  width = 8, height = 2)


featcor = read.csv(paste0(datadir, 'featcor_toplot.csv'), row.names = 1)
paletteLength <- 200
heat_scale <- colorRampPalette(c("#003836", "white", "#a87b00"))(paletteLength)

myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1,
                  length.out=floor(paletteLength/2)))
pdf(paste0(figdir, "corheatmap-allfeats.pdf"),
    width=6.5, height=6)
ph = pheatmap(featcor, color = heat_scale,
         breaks = myBreaks, show_colnames = F,
         show_rownames = F)
dev.off()
#save_pheatmap_pdf(ph, filename = paste0(figdir, "corheatmap-allfeats.pdf"))

# # just schematics for compound descriptor generation
# plate1= read.csv(paste0(datadir, 'venetoclax-180528_Plate5.csv'), row.names = 1)
# col_clust = pheatmap(cor(plate1, use="pairwise.complete.obs"), silent = T)
# plate1_annot = read.csv(paste0(datadir, 'annot-venetoclax-180528_Plate5.csv'), row.names = 1)
# plate1 = plate1[which(plate1_annot$Culture == 'Mono-culture'),]
# 
# row_clust = pheatmap(cor(t(plate1), use="pairwise.complete.obs"), silent = T)
# 
# ph = pheatmap(plate1,
#          breaks = myBreaks,
#          cluster_rows = row_clust$tree_row, 
#          cluster_cols = col_clust$tree_col,
#          color=heat_scale,
#          show_colnames = F, show_rownames = F,
#          treeheight_row = 0, 
#          treeheight_col = 0)
# save_pheatmap_pdf(ph, filename = paste0(figdir, 'schematic1.pdf'),
#                   width = 3, height = 4)
# 
# medprof = t(as.matrix(matrixStats::colMedians(as.matrix(plate1))))
# colnames(medprof) = colnames(plate1)
# ph = pheatmap(medprof,
#          breaks = myBreaks,
#          cluster_rows = F, 
#          cluster_cols = col_clust$tree_col,
#          color=heat_scale,
#          show_colnames = F, show_rownames = F,
#          treeheight_row = 0, 
#          treeheight_col = 0)
# save_pheatmap_pdf(ph, filename = paste0(figdir, 'schematic1.pdf'),
#                   width = 3, height = 0.2)
# 
# medprof = (as.matrix(drugprof_sel['Venetoclax_0.225',]))
# colnames(medprof) = colnames(drugprof_sel)
# ph = pheatmap(medprof,
#               breaks = myBreaks,
#               cluster_rows = F, 
#               cluster_cols = col_clust$tree_col,
#               color=heat_scale,
#               show_colnames = F, show_rownames = F,
#               treeheight_row = 0, 
#               treeheight_col = 0)

