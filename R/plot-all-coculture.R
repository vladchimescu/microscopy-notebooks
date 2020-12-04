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
set.seed(512)
drugs_exclude = c("SSZ", "PEITC",
                  "SSZ + doxorubicine",
                  "SSZ + fludarabine",
                  "Vorinostat + PEITC",
                  "Hydroxychloroquine",
                  "β-phenylethyl isothiocyanate")

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

drug_prof = drug_prof[!rownames(drug_prof) %in% drugs_exclude,]

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

drugs_autophagy = c("Bafilomycin A1",
                    "Ganetespib",
                    "JQ1",
                    "UMI-77",
                    "Venetoclax",
                    "Cytarabine",
                    "Fludarabine",
                    "Ixazomib",
                    "Vorinostat",
                    "Carfilzomib")
df = read.delim(paste0(datadir, 'coculture_selected_feats_all.txt'))
df = filter(df, !Drug %in% drugs_exclude)
feats = c("ch-Lysosomal-area", 
          "ch-Lysosomal-mean_intensity")
lys_df = filter(df, feature %in% feats) %>%
  filter(Drug %in% drugs_autophagy)
# lys_df = filter(inner_join(lys_df, patannot),
#                 !(feature == 'ch-Lysosomal-mean_intensity' & Culture == 'Co-culture' & Diagnosis == 'CLL'))
library(ggplot2)
# ggplot(lys_df,
#        aes(x = Drug, y= value, fill = Culture))+
#   geom_boxplot()+
#   facet_wrap(~ feature, ncol=2) +
#   scale_fill_manual(values = c("#407088", "#b3b3b3")) + 
#   theme_bw() + 
#   theme(aspect.ratio = 1)

p = ggplot(filter(lys_df, feature == 'ch-Lysosomal-area'),
       aes(x = Drug, y= value, fill = Culture))+
  geom_boxplot()+
  scale_fill_manual(values = c("#407088", "#b3b3b3")) + 
  labs(x = "", y = "Normalized lysosomal area") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = 'none',
        legend.text = element_text(size=14),
        legend.title = element_text(size=16)) 
ggsave(p, filename = paste0(figdir, 'drugs-lysosomal-area.pdf'),
       width = 6, height = 5)

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



# lysosomal area by diagnosis
feat_df = inner_join(df, viab_df)
feat_df = filter(feat_df, Drug %in% drugs_autophagy& feature %in% c("ch-Lysosomal-area")) %>%
  mutate(feature = plyr::mapvalues(feature, from = c("ch-Lysosomal-area"),
                                   to = c("Lysosomal area")))
feat_df = inner_join(feat_df, patannot) %>%
  mutate(Diagnosis = ifelse(Diagnosis == 'CLL', 'CLL samples', 'non-CLL entities'))
p = ggplot(feat_df,
           aes(x = n_viab_norm, y = value, color = Culture))+
  geom_point() +
  scale_color_manual(values = c("#407088", "#b3b3b3")) + 
  ggpubr::stat_cor(aes(label = ..r.label..),
                   label.y.npc= 1, 
                   label.x.npc = 0.6,
                   method = "pearson",
                   cor.coef.name = 'R', 
                   show.legend = F) + 
  facet_wrap(~ Diagnosis + Culture, ncol= 2, scales = 'free') +
  labs(x = "Normalized viable cell count",
       y = "Normalized lysosomal area") + 
  theme_classic() +
  #coord_fixed() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        legend.position = 'none')
ggsave(p, filename = paste0(figdir, 'lysosomal-area-vs-cell-count.png'),
       width = 6, height = 6)


feat_df = inner_join(df, viab_df)
feat_df = filter(feat_df, Drug %in% drugs_autophagy& feature %in% c("ch-Lysosomal-mean_intensity")) %>%
  mutate(feature = plyr::mapvalues(feature, from = c("ch-Lysosomal-mean_intensity"),
                                   to = c("Lysosomal intensity")))

feat_df = inner_join(feat_df, patannot) %>%
  mutate(Diagnosis = ifelse(Diagnosis == 'CLL', 'CLL samples', 'non-CLL entities'))
p = ggplot(feat_df,
           aes(x = n_viab_norm, y = value, color = Culture))+
  geom_point() +
  scale_color_manual(values = c("#407088", "#b3b3b3")) + 
  ggpubr::stat_cor(aes(label = ..r.label..),
                   label.y.npc= 1, 
                   label.x.npc = 0.7,
                   method = "pearson",
                   cor.coef.name = 'R', 
                   show.legend = F) + 
  facet_wrap(~ Diagnosis + Culture, ncol= 2, scales = 'free') +
  labs(x = "Normalized viable cell count",
       y = "Normalized lysosomal intensity") + 
  theme_classic() +
  #coord_fixed() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        legend.position = 'none')
ggsave(p, filename = paste0(figdir, 'lysosomal-intensity-vs-cell-count.png'),
       width = 7, height = 7)

rep1 = filter(df, plate == '181109_Plate1')
rep2 = filter(df, plate == '181109_Plate4')

cor_df = inner_join(rep1, rep2,
           by = c("well", "Drug", 
                  "conc", "Culture", "feature"))

highcor = group_by(cor_df, feature) %>%
  summarise(R = cor(value.x, value.y)) %>%
  filter(R > 0.5)

cor_df = filter(cor_df, feature %in% unique(highcor$feature))

pdf(paste0(figdir,"nonCLL-featcor-by-coculture.pdf"), width=16, height=16)
pages = ceiling(length(unique(cor_df$feature))/16)
for (i in seq_len(pages)) {
  page = ggplot(cor_df,
         aes(x=value.x, y=value.y, color = Culture))+
    geom_point(alpha=0.5) +
    scale_color_manual(values = c("#407088", "#b3b3b3")) + 
    ggpubr::stat_cor(cor.coef.name = 'R', show.legend = F) + 
    facet_wrap_paginate(~ feature, ncol= 4, nrow = 4, page=i, scales = "free") +
    labs(x = "biol rep 1",
         y = "biol rep 2") + 
    theme_bw() + 
    theme(aspect.ratio = 1)
  print(page)   
}
dev.off()


p = ggplot(filter(cor_df, feature  %in% feats) %>%
         mutate(feature = plyr::mapvalues(feature,
                                          from=feats,
                                          to=c("Lysosomal area",
                                               "Lysosomal intensity"))),
       aes(x=value.x, y=value.y))+
  geom_point(alpha=0.3, color="#ec5858") +
  facet_wrap(~ feature + Culture, ncol = 2) +
  geom_abline(slope=1, linetype='dotted') + 
  ggpubr::stat_cor() +
  labs(x = "AML (biol rep 1)",
       y = "AML (biol rep 2)") +
  theme_bw() +
  theme(aspect.ratio = 1)
ggsave(p, filename = paste0(figdir, 'AML-biolrep-lysosomal-by-culture.pdf'),
       width = 6, height = 6)

# now do this for CLL replicates
patannot = readxl::read_xlsx("~/Documents/embl/gitlab/microscopy/data/Metafiles/Screened_Patients.xlsx")
patannot = mutate(patannot, Screen = gsub("181106", "181109", Screen),
                  Plate_Nr = gsub("_", "", Plate_Nr)) %>%
  mutate(plate = paste(Screen, Plate_Nr, sep = "_")) %>%
 select(plate, OMZ_Pat_ID)

rep1 = bind_rows(filter(df, plate == '180306_Plate1'),
                 filter(df, plate == '180306_Plate3'),
                 filter(df, plate == '180504_Plate6')) %>%
  inner_join(patannot)
rep2 = bind_rows(filter(df, plate == '180306_Plate5'),
                 filter(df, plate == '180424_Plate6'),
                 filter(df, plate == '180504_Plate10')) %>%
  inner_join(patannot)

cor_df = inner_join(rep1, rep2,
                    by = c("OMZ_Pat_ID","well", "Drug", 
                           "conc", "Culture", "feature"))

p = ggplot(filter(cor_df, feature  %in% feats) %>%
             mutate(feature = plyr::mapvalues(feature,
                                              from=feats,
                                              to=c("Lysosomal area",
                                                   "Lysosomal intensity"))),
           aes(x=value.x, y=value.y))+
  geom_point(alpha=0.3, color="#ec5858") +
  facet_wrap(~ feature + Culture, ncol = 2) +
  geom_abline(slope=1, linetype='dotted') + 
  ggpubr::stat_cor() +
  labs(x = "CLL (biol rep 1)",
       y = "CLL (biol rep 2)") +
  theme_bw() +
  theme(aspect.ratio = 1)
ggsave(p, filename = paste0(figdir, 'CLL-biolrep-lysosomal-by-culture.pdf'),
       width = 6, height = 6)


highcor = group_by(cor_df, feature) %>%
  summarise(R = cor(value.x, value.y)) %>%
  filter(R > 0.5)

cor_df = filter(cor_df, feature %in% unique(highcor$feature))
pdf(paste0(figdir,"CLL-featcor-by-coculture.pdf"), width=16, height=16)
pages = ceiling(length(unique(cor_df$feature))/16)
for (i in seq_len(pages)) {
  page = ggplot(cor_df,
                aes(x=value.x, y=value.y, color = Culture))+
    geom_point(alpha=0.5) +
    scale_color_manual(values = c("#407088", "#b3b3b3")) + 
    ggpubr::stat_cor(cor.coef.name = 'R', show.legend = F) + 
    facet_wrap_paginate(~ feature, ncol= 4, nrow = 4, page=i, scales = "free") +
    labs(x = "biol rep 1",
         y = "biol rep 2") + 
    theme_bw() + 
    theme(aspect.ratio = 1)
  print(page)   
}
dev.off()

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

