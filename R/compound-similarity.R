## Compound similarity based on profiles
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggnewscale)

datadir = "~/Documents/github/microscopy-notebooks/data/"
figdir = "~/Documents/github/microscopy-notebooks/figures/compclust/"
if(!dir.exists(figdir)) dir.create(figdir)

# global variable: number of principal components
npc <- 20

cormat <- function(mat, thresh) {
  #' Correlation matrix based on first PCs
  #' 
  #' Computes the correlation matrix based on the first npc
  #' principal components and returns only correlations above
  #' the specified threshold
  #' 
  #' @param mat Data matrix
  #' @param thresh Correlation threshold
  #' 
  #' @return Thresholded correlation matrix
  pr = prcomp(mat)
  cormat = cor(t(pr$x[,1:npc]))
  cormat[cormat < thresh] = 0
  cormat
}

extract_triangular <- function(prof_mono, prof_co, levs, thresh = 0.5) {
  #' Extract trinagular matrices from mono- and coculture profile similarities
  #' 
  #' Combines the upper triangular matrix of coculture profile similarity
  #' with the lower triangular matrix of monoculture profile similarity
  #' 
  #' @param prof_mono Monoculture profiles
  #' @param prof_co Coculture profiles
  #' @param levs Drug order
  #' @param thresh Correlation threshold. Passed to @method cormat
  #' 
  #' @return List with monoculture and coculture profile similarities
  cormat_mono = cormat(prof_mono, thresh = thresh)
  cormat_co = cormat(prof_co, thresh = thresh)
  
  cormat_co = cormat_co[levs, levs]
  cormat_mono = cormat_mono[levs, levs]
  
  cormat_co[upper.tri(cormat_co, diag=F)] = NA
  cormat_mono[lower.tri(cormat_mono, diag=T)] = NA
  
  heat_co = reshape2::melt(cormat_co, varnames = c("row", "col"),
                           value.name = "Coculture") %>%
    filter(!is.na(Coculture))
  heat_mono =  reshape2::melt(cormat_mono, varnames = c("row", "col"),
                              value.name = "Monoculture") %>%
    filter(!is.na(Monoculture))
  
  list(cocult = heat_co,
       monocult = heat_mono)
}

plot_similarity <- function(heat_list, lims = c(0,1), textsize = 9) {
  #' Plot profile similarity
  #' 
  #' Heatmap of coculture and monoculture profile similarity
  #' 
  #' @param heat_list List with mono- and coculture profile similarities
  #' @param lims Fill scale limits (optional)
  #' @param textsize Base font size (optional)
  #' 
  #' @return ggplot object
  p = ggplot(heat_list$cocult, 
             aes(x = col,
                 y = row, 
                 fill = Coculture)) +
    geom_tile(color = "white") +
    labs(x = "", y = "") + 
    scale_fill_gradient(high = "#ee854a",
                        low = "white",
                        limits=c(0,1)) + 
    new_scale("fill") + 
    geom_tile(data = heat_list$monocult, 
              aes(x = col,
                  y = row,  
                  fill = Monoculture),
              inherit.aes = F,
              color = "white") +
    scale_fill_gradient(high = "#4878d0",
                        low = "white",
                        limits=lims) +
    scale_y_discrete(position = 'right') + 
    theme_classic(base_size = textsize) +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5,
                                     hjust = 1),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 1.25),
          axis.line = element_blank(),
          axis.ticks = element_blank())+
    coord_fixed(ratio = 1)
  
  p
}

set.seed(274)
drugs_exclude = c("SSZ", "PEITC",
                  "SSZ + doxorubicine",
                  "SSZ + fludarabine",
                  "Vorinostat + PEITC",
                  "Hydroxychloroquine",
                  "Oxaliplatin",
                  "β-phenylethyl isothiocyanate")

cll_prof = readRDS(paste0(datadir, 'CLL-profiles-wide.rds'))
noncll_prof = readRDS(paste0(datadir, 'nonCLL-profiles-wide.rds'))
stopifnot(all(rownames(cll_prof) == rownames(noncll_prof)))
drug_prof = cbind(cll_prof, noncll_prof)

drug_prof = drug_prof[!rownames(drug_prof) %in% drugs_exclude,]
drug_prof = drug_prof[!grepl('Vitamin D', rownames(drug_prof)),]

rownames(drug_prof) = plyr::mapvalues(rownames(drug_prof),
                                      from = c("IRAK1/4 Inhibitor I",
                                               "IRAK4 Inhibitor, Compound 26"),
                                      to = c("IRAK1 inhibitor",
                                             "IRAK4 inhibitor"))
prof_mono = drug_prof[,grepl("_M", colnames(drug_prof))]
prof_co = drug_prof[,grepl("_C", colnames(drug_prof))]

# cluster by all profiles and compute similairty matrices
# for monoculture and coculture separately
row_clust = pheatmap(drug_prof, silent = T)
levs = row_clust$tree_row$labels[row_clust$tree_row$order]
heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs,
                               thresh = 0.4)
p = plot_similarity(heat_list)
ggsave(p + theme(legend.position = "none"),
       filename = paste0(figdir, "compound-similarity-", npc, "-pcs-thresh-clust-by-all-profiles.pdf"),
       width = 13.5, height = 13.5, units = 'cm')

# now try to cluster by dimension-reduced compound profiles
pr = prcomp(drug_prof)
row_clust = pheatmap(pr$x[,1:npc], silent = T)
levs = row_clust$tree_row$labels[row_clust$tree_row$order]

pdf(paste0(figdir, "hclust-drugs.pdf"), height = 6, width = 2)
plot(as.dendrogram(row_clust$tree_row), leaflab='none', horiz = T, xaxt = 'n', yaxt = 'n')
dev.off()

heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs,
                               thres = 0.4)
p = plot_similarity(heat_list)
ggsave(p + theme(legend.position = "none"),
       filename = paste0(figdir, "compound-similarity-", npc, "-pcs-thresh-clust-by-all-profiles-dimred.pdf"),
       width = 13.5, height = 13.5, units = 'cm')

# now remove the cell count from compound profiles => only image profiles
drug_prof = drug_prof[,!grepl("cellcount", colnames(drug_prof))]
prof_mono = drug_prof[,grepl("_M", colnames(drug_prof))]
prof_co = drug_prof[,grepl("_C", colnames(drug_prof))]

row_clust = pheatmap(drug_prof, silent = T)
levs = row_clust$tree_row$labels[row_clust$tree_row$order]
heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs,
                               thresh = 0.4)
p = plot_similarity(heat_list)
ggsave(p + theme(legend.position = "none"),
       filename = paste0(figdir, "img-profile-similarity-", npc, "-pcs-thresh-clust-by-image-profiles.pdf"),
       width = 13.5, height = 13.5, units = 'cm')


# now try to cluster by dimension-reduced image profiles
pr = prcomp(drug_prof)
row_clust = pheatmap(pr$x[,1:npc], silent = T)
levs = row_clust$tree_row$labels[row_clust$tree_row$order]

heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs, thresh = 0.4)
p = plot_similarity(heat_list)
ggsave(p + theme(legend.position = "none"),
       filename = paste0(figdir, "img-profile-similarity-", npc, "-pcs-thresh-clust-by-image-profiles-dimred.pdf"),
       width = 13.5, height = 13.5, units = 'cm')

# save the legend
leg = as_ggplot(get_legend(p + theme(legend.position = 'top',
                                             legend.key.width = unit(5, 'mm'),
                                             legend.key.height = unit(2,'mm'))))
ggsave(leg, filename = paste0(figdir, 'similarity-legend.pdf'),
       width=10, height = 2, units = 'cm')

# now cluster by coculture profiles
pr = prcomp(prof_co)
row_clust = pheatmap(pr$x[,1:npc], silent = T)
levs = row_clust$tree_row$labels[row_clust$tree_row$order]
heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs,
                               thresh = 0.4)
p = plot_similarity(heat_list)
ggsave(p + theme(legend.position = "none"),
       filename = paste0(figdir, "compound-similarity-", npc, "-pcs-thresh-clust-by-coculture-profiles.pdf"),
       width = 13.5, height = 13.5, units = 'cm')


# now cluster by monoculture profiles
pr = prcomp(prof_mono)
row_clust = pheatmap(pr$x[,1:npc], silent = T)
levs = row_clust$tree_row$labels[row_clust$tree_row$order]
heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs,
                               thresh = 0.4)
p = plot_similarity(heat_list)
ggsave(p + theme(legend.position = "none"),
       filename = paste0(figdir, "compound-similarity-", npc, "-pcs-thresh-clust-by-monoculture-profiles.pdf"),
       width = 13.5, height = 13.5, units = 'cm')


heat_list = extract_triangular(prof_mono,
                               prof_co,
                               levs,
                               thres = -10)

cocult_sim = filter(heat_list$cocult, row != col) %>%
  mutate(drugcomb = purrr::map2_chr(as.character(row),
                                    as.character(col),
                                    ~ paste(sort(c(.x, .y)), 
                                            collapse = "_"))) %>%
  select(drugcomb, Coculture) %>%
  rename(cor = Coculture)

monocult_sim = filter(heat_list$monocult, row != col) %>%
  mutate(drugcomb = purrr::map2_chr(as.character(row),
                                    as.character(col),
                                    ~ paste(sort(c(.x, .y)), 
                                            collapse = "_"))) %>%
  select(drugcomb, Monoculture) %>%
  rename(cor = Monoculture)

tanimoto = read.csv(paste0(datadir,
                               "coculture_metafiles/FCPF4-jaccard-dist.csv"),
                    row.names = 1,
                        check.names = F)
tanimoto_df = reshape2::melt(as.matrix(tanimoto),
                             value.name = "tanimoto",
                             varnames = c("d1", "d2"))
tanimoto_df = mutate(tanimoto_df, d1 = as.character(d1),
                     d2 = as.character(d2))
tanimoto_df = filter(tanimoto_df, d1 != d2)

tanimoto_df = mutate(tanimoto_df, drugcomb = purrr::map2_chr(as.character(d1),
                                                as.character(d2),
                                                ~ paste(sort(c(.x, .y)), 
                                                        collapse = "_"))) %>%
  distinct(drugcomb, tanimoto) %>%
  mutate(sim = 1 - tanimoto)

sim_df = inner_join(monocult_sim, cocult_sim, by='drugcomb') %>%
  inner_join(tanimoto_df)
p = ggplot(filter(sim_df, sim < 0.5),
       aes(x = cor.x, y = cor.y, color = sim))+
  geom_point(alpha = 0.75) +
  geom_point(data = filter(sim_df, sim >= 0.5 ), alpha = 0.75) + 
  #scale_color_viridis_c(option = 'A') + 
  scale_colour_gradientn(colours = c("white", "orange", "red"),
                         values = c(0, 0.3, 0.6, 1)) + 
  labs(x = "Profile similarity in monoculture",
       y = "Profile similarity in coculture",
       color = "Structural similarity") + 
  geom_text(data = filter(sim_df, sim > 0.5 & (cor.x > 0.5 | cor.y > 0.5 | cor.x < 0) ), 
            aes(label = drugcomb), size = 7 * 5/ 14) + 
  theme_classic(base_size = 9) +
  stat_cor(label.y = 1)
ggsave(p + theme(legend.position = 'none'), 
       filename = paste0(figdir, "coculture-monoculture-sim-scatter.pdf"),
       width = 7.5, height = 7.5, units = 'cm')

# save the legend
leg = as_ggplot(get_legend(p ))
ggsave(leg, filename = paste0(figdir, 'profile-chemsimilarity-legend.pdf'),
       width=3, height = 5, units = 'cm')
