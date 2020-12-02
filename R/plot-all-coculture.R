## plot all coculture profiles together
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


cll_prof = readRDS(paste0(datadir, 'CLL-profiles-wide.rds'))
noncll_prof = readRDS(paste0(datadir, 'nonCLL-profiles-wide.rds'))

stopifnot(all(rownames(cll_prof) == rownames(noncll_prof)))

drug_prof = cbind(cll_prof, noncll_prof)
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
  select(plate, Diagnosis)
colannot = left_join(colannot, patannot)

rownames(colannot) = colannot$names
colannot = select(colannot, -c(names, plate))

chan = c("#63b7af", "#035aa6", "#d8345f")
names(chan) = c("Calcein", "Hoechst", "Lysosomal")
dg = c("#010059", "#107595", "#fdbfb3", "#fcf594", "#aee7e8")
names(dg) = c("AML", "CLL", "HCL", "MCL", "T-PLL")
# for culture
colscult = c("#407088", "#b3b3b3")
names(colscult) = c("Coculture", "Monoculture")
ann_colors <- list(Channel=chan,
                   Diagnosis=dg,
                   Culture=colscult)

rownames(drug_prof) = plyr::mapvalues(rownames(drug_prof),
                                      from = c("β-phenylethyl isothiocyanate",
                                               "IRAK4 Inhibitor, Compound 26"),
                                      to = c("PEITC",
                                             "IRAK4 Inhibitor"))

paletteLength <- 200
heat_scale <- colorRampPalette(c("#003836", "white", "#a87b00"))(paletteLength)
maxval = 6
myBreaks <- c(seq(-maxval, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, maxval,
                  length.out=floor(paletteLength/2)))
row_clust = pheatmap(cor(t(drug_prof), use="pairwise.complete.obs"), silent = T)
#col_clust = pheatmap(cor(drug_prof, use="pairwise.complete.obs"), silent = T)

# plot the heatmap with all 57 conditions
ph = pheatmap(drug_prof,
              cluster_rows = row_clust$tree_row, 
              #cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 
                                        'drugprofiles-coculture-all-aggregated.pdf'),
                  width = 9, height = 8)


df = read.delim(paste0(datadir, 'coculture_selected_feats_all.txt'))
feats = c("ch-Lysosomal-area", 
          "ch-Lysosomal-mean_intensity")
lys_df = filter(df, feature %in% feats) %>%
  filter(Drug %in% c("Bafilomycin A1",
                     "Ganetespib",
                     "Venetoclax",
                     "JQ1",
                     "Vorinostat",
                     "Carfilzomib",
                     "Hydroxychloroquine"))
library(ggplot2)
ggplot(lys_df, aes(x = Drug, y= value, fill = Culture))+
  geom_boxplot()+
  facet_wrap(~ feature, ncol=2) +
  theme_bw() + 
  theme(aspect.ratio = 1)

normalize_viab <- function(viab_df) {
  ctrl_df = filter(viab_df, Drug== 'DMSO')
  ctrl_df = group_by(ctrl_df, plate, Culture) %>%
    # remove control wells at the edge
    filter(Edge < 2) %>%
    summarise(mu = mean(n_viab),
              sd = sd(n_viab))
  
  viab_df = inner_join(viab_df, ctrl_df)
  viab_df = mutate(viab_df, n_viab_norm = (n_viab - mu) / sd)
  viab_df
}

select_cols <- function(viab_df) {
  select(viab_df, plate, well, n_viab_norm)
}

# load viability data
cll_viab = readRDS('~/Documents/embl/gitlab/microscopy/data/CLL_viab_2020.rds')
cll_viab = normalize_viab(cll_viab)
noncll_viab = readRDS('~/Documents/embl/gitlab/microscopy/data/AML_viab_2020.rds') %>%
  rename(n_viab = viable)
noncll_viab = normalize_viab(noncll_viab)

cll_viab = select_cols(cll_viab)
noncll_viab = select_cols(noncll_viab)

viab_df = bind_rows(cll_viab, noncll_viab)

# intersect with viability data
feat_df = inner_join(df, viab_df)

sums = rowSums(abs(drug_prof), na.rm = T)
names(sums) = rownames(drug_prof)
top20 = names(sort(sums, decreasing = T)[1:20])
feat_df = filter(feat_df, Drug %in% top20)

library(ggforce)
pdf(paste0(figdir,"cor-viab-feature-coculture.pdf"), width=16, height=16)
pages = ceiling(length(unique(feat_df$feature))/16)
for (i in seq_len(pages)) {
  page <- ggplot(feat_df,
                 aes(x = n_viab_norm, y = value, color = Culture))+
    geom_point(alpha = 0.3) +
    scale_color_manual(values = c("#407088", "#b3b3b3")) + 
    ggpubr::stat_cor(cor.coef.name = 'R', show.legend = F) + 
    facet_wrap_paginate(~ feature, ncol= 4, nrow = 4, page=i, scales = "free") +
    theme_bw() +
    theme(aspect.ratio = 1)
  print(page)   
}
dev.off()


feat_df = inner_join(df, viab_df) %>%
  filter(Drug != 'DMSO')
sums = rowSums(abs(drug_prof), na.rm = T)
names(sums) = rownames(drug_prof)
# top20 = names(sort(sums, decreasing = T)[1:20])
# feat_df = filter(feat_df, Drug %in% top20)
p = ggplot(filter(feat_df, feature %in% feats) %>%
         mutate(feature = plyr::mapvalues(feature, from = feats,
                                          to = c("Lysosomal area",
                                                 "Lysosomal intensity"))),
       aes(x = n_viab_norm, y = value, color = Culture))+
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("#407088", "#b3b3b3")) + 
  ggpubr::stat_cor(aes(label = ..r.label..),
                   method = "pearson",cor.coef.name = 'R', show.legend = F) + 
  #ggpubr::stat_cor(method = "spearman",cor.coef.name = 'rho', show.legend = F) + 
  facet_wrap(~ feature, ncol= 4, nrow = 4) +
  labs(x = "Standardized viable cell count",
       y = "Standardized value") + 
  theme_bw() +
  coord_fixed() +
  theme(aspect.ratio = 1,
        legend.position = 'none')
ggsave(p, filename = paste0(figdir, 'lysosomal-vs-cell-count.png'),
       width = 6, height = 3.5)

# # for spontaneous apoptosis center and scale co- and monoculture wells jointly
# normalize_control_viab <- function(viab_df) {
#   ctrl_df = filter(viab_df, Drug== 'DMSO')
#   ctrl_df = group_by(ctrl_df, plate) %>%
#     # remove control wells at the edge
#     filter(Edge < 2) %>%
#     summarise(mu = mean(n_viab),
#               sd = sd(n_viab))
#   
#   viab_df = inner_join(viab_df, ctrl_df)
#   viab_df = mutate(viab_df, n_viab_norm = (n_viab - mu) / sd)
#   filter(viab_df, Drug == 'DMSO')
# }
# 
# cll_viab = readRDS('~/Documents/embl/gitlab/microscopy/data/CLL_viab_2020.rds')
# cll_viab = normalize_control_viab(cll_viab)
# noncll_viab = readRDS('~/Documents/embl/gitlab/microscopy/data/AML_viab_2020.rds') %>%
#   rename(n_viab = viable)
# noncll_viab = normalize_control_viab(noncll_viab)
# 
# cll_viab = select_cols(cll_viab)
# noncll_viab = select_cols(noncll_viab)
# ctrl_df = bind_rows(cll_viab, noncll_viab)
# 
# # explore control wells
# feat_df = inner_join(df, ctrl_df) %>%
#   filter(Drug == 'DMSO')
# 
# pdf(paste0(figdir,"cor-spapoptosis-feature-coculture.pdf"), width=16, height=16)
# pages = ceiling(length(unique(feat_df$feature))/16)
# for (i in seq_len(pages)) {
#   page <- ggplot(feat_df,
#                  aes(x = n_viab_norm, y = value, color = Culture))+
#     geom_point(alpha = 0.3) +
#     scale_color_manual(values = c("#407088", "#b3b3b3")) + 
#     ggpubr::stat_cor(cor.coef.name = 'R', show.legend = F) + 
#     facet_wrap_paginate(~ feature, ncol= 4, nrow = 4, page=i, scales = "free") +
#     theme_bw() +
#     theme(aspect.ratio = 1)
#   print(page)   
# }
# dev.off()

