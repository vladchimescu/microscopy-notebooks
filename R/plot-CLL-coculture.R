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

drug_prof = read.csv(paste0(datadir, 'all_CLL_profiles.csv'),
                     row.names = 1)

colannot = data.frame(names = colnames(drug_prof))
colannot = mutate(colannot, Channel = ifelse(grepl('Calcein', names), 'Calcein',
                                             ifelse(grepl('Hoechst', names), 'Hoechst',
                                                    ifelse(grepl('Lysosomal', names),
                                                           'Lysosomal', 'Cell count'))))
colannot = mutate(colannot, Culture = ifelse(grepl("_C\\.", names), "Coculture", "Monoculture"))
colannot = mutate(colannot, plate = gsub("(.+_[CM]\\.)(18.+)", '\\2', names))

# add diagnosis based on the screened plate
patannot = readxl::read_xlsx("~/Documents/embl/gitlab/microscopy/data/Metafiles/Screened_Patients.xlsx")
patannot = mutate(patannot, Screen = gsub("181106", "181109", Screen),
         Plate_Nr = gsub("_", "", Plate_Nr)) %>%
  mutate(plate = paste(Screen, Plate_Nr, sep = "_")) %>%
  select(plate, IGHV)
colannot = left_join(colannot, patannot)
colannot = mutate(colannot, IGHV = ifelse(is.na(IGHV), 'NA', IGHV))

rownames(colannot) = colannot$names
colannot = select(colannot, -c(names, plate))

chan = c("#035aa6", "#c20078", "#ffd66b")
names(chan) = c("Hoechst", "Lysosomal", "Cell count")
dg = c("#1e488f", "#d46a7e", "#ececec")
names(dg) = c('U', 'M', 'NA')
# for culture
colscult = c('#ee854a', '#4878d0')
names(colscult) = c("Coculture", "Monoculture")
ann_colors <- list(Channel=chan,
                   IGHV=dg,
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

save_pheatmap_pdf(ph, filename = paste0(figdir, 'CLL-drugprofiles-all.pdf'),
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
save_pheatmap_pdf(ph, filename = paste0(figdir, 'CLL-drugprofiles-selected.pdf'),
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
save_pheatmap_pdf(ph, filename = paste0(figdir, 'CLL-drugprofiles-small-subset.pdf'),
                  width = 9, height = 8)

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
save_pheatmap_pdf(ph, filename = paste0(figdir, 'CLL-drugprofiles-autophagy.pdf'),
                  width = 9, height = 8)

# Coculture vs Monoculture morphological feature difference by sample
cultdiff = read.csv(paste0(datadir, 'CLL_diff_DMSO.csv'))
cultdiff = dplyr::rename(cultdiff, feature=X)
cultdiff = filter(cultdiff, feature %in% c("ch-Hoechst-eccentricity",
                                           "ch-Hoechst-InfoMeas1-d5-3",
                                           "ch-Lysosomal-area",
                                           "ch-Lysosomal-extent")) %>%
  mutate(feature = plyr::mapvalues(feature, from = c("ch-Hoechst-eccentricity",
                                                     "ch-Hoechst-InfoMeas1-d5-3",
                                                     "ch-Lysosomal-area",
                                                     "ch-Lysosomal-extent"),
                                   to = c("Hoechst eccentricity",
                                          "Hoechst InfoMeas1 [d = 5]",
                                          "Lysosomal area",
                                          "Lysosomal extent")))

dmsodiff_out = mutate(cultdiff, Diagnosis = 'CLL') %>%
  group_by(Diagnosis, feature) %>%
  summarise(mean_diff_medians = mean(diff_medians)) %>%
  ungroup()

saveRDS(dmsodiff_out, file = paste0(datadir, "CLL-mean_diff_medians.rds"))

df_wide = tidyr::spread(select(cultdiff, feature, plate, diff_medians),
                        key='feature', value='diff_medians')
rownames(df_wide) = df_wide$plate
df_wide = df_wide[,-1]
colnames(df_wide) = gsub("ch-", "", colnames(df_wide))

df_wide = t(df_wide)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#417ca8","white", "#d1483a"))(paletteLength)
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
ann_colors <- list(IGHV=dg)

df_wide[df_wide == 0] = NA

ph = pheatmap(df_wide,
              cluster_rows = F,
              cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              fontsize = 9,
              treeheight_row = 0, 
              treeheight_col = 0,
              na_col = '#888888')

save_pheatmap_pdf(ph, filename = paste0(figdir, 'CLL-DMSO-culture-difference.pdf'),
                  width = 17.5/2.54, height = 2.5/2.54)

# aggregate the profiles across concentrations
drug_prof$drugconc = rownames(drug_prof)
drugprof_long = reshape2::melt(drug_prof, id.variable='drugconc') %>%
  tidyr::separate(drugconc, c('drug', 'conc'), sep='_') %>%
  filter(drug != "DMSO")

# take the extreme normalized value for each feature across 3 concentrations
drugprof_long = group_by(drugprof_long, drug, variable) %>%
  summarise(value = value[which.max(abs(value))])

heat_wide = tidyr::spread(drugprof_long, variable, value)
mat_heatmap = as.matrix(heat_wide[,-1])
rownames(mat_heatmap) = heat_wide$drug
saveRDS(mat_heatmap, file = paste0(datadir, 'CLL-profiles-wide.rds'))
