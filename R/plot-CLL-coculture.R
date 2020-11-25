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
                                             ifelse(grepl('Hoechst', names), 'Hoechst', 'Lysosomal')))
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

chan = c("#035aa6", "#d8345f")
names(chan) = c("Hoechst", "Lysosomal")
dg = c("#85a392", "#f2a365", "#ececec")
names(dg) = c('U', 'M', 'NA')
# for culture
colscult = c("#407088", "#b3b3b3")
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
