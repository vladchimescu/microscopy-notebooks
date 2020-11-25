library(RSQLite)
library(DBI)
library(dplyr)
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

#----------Explore single cell profiles--------------
con <- dbConnect(SQLite(), paste0(datadir,
                                  "Bray-scprofiles/24277.sqlite"))
as.data.frame(dbListTables(con))

cellfeat <- dbGetQuery(con, "SELECT * FROM Cells LIMIT 1000")
image = dbReadTable(con, 'Image')

dbGetQuery(con, 'SELECT Count(*) FROM Image')
dbGetQuery(con, 'SELECT Count(*) FROM Cells')
# explore feature names
feats = colnames(cellfeat)

#----------Load Metadata and plate annotation-------------------------------------------
chemannot = read.csv(paste0(datadir,'Bray-metadata/chemical_annotations.csv'))
barcodes = read.csv(paste0(datadir,'Bray-metadata/barcode_platemap.csv'))

maps = grep(".txt", list.files(paste0(datadir, 'Bray-metadata/platemap')), value = T)

platemap = lapply(maps, function(x) read.table(paste0(datadir, 'Bray-metadata/platemap/', x), sep = '\t',
                                               header = T))
names(platemap) = gsub("\\.txt", "", maps)
platemap = plyr::ldply(platemap, .id='Plate_Map_Name')

platemap = inner_join(barcodes, platemap) %>%
  rename(plateID = Assay_Plate_Barcode)

platemap = select(platemap, Plate_Map_Name, plateID, 
       well_position, ASSAY_WELL_ROLE,
       broad_sample, mmoles_per_liter)

chemannot = select(chemannot, BROAD_ID, CPD_NAME,
       CPD_NAME_TYPE, CPD_SMILES)
chemannot = rename(chemannot, broad_sample = BROAD_ID)
platemap = left_join(platemap, chemannot)

# all drugs
alldrugs = distinct(platemap, CPD_NAME, CPD_NAME_TYPE)

write.table(platemap, file = paste0(datadir, 'Bray-metadata/plate_annot.txt'),
            sep = '\t',
            col.names = T,
          row.names = F, quote = F)


# only treated wells, check out the replicates
drugs = filter(platemap, ASSAY_WELL_ROLE == 'treated')


#---------Mean well profiles-----------
drugprof = read.csv(paste0(datadir, 'drugprofiles_Bray.csv'))
rownames(drugprof) = paste(drugprof$CPD_NAME, drugprof$mmoles_per_liter, sep = '_')
drugprof = select(drugprof, -c(CPD_NAME, mmoles_per_liter))

colannot = data.frame(Level = as.character(sapply(colnames(drugprof),
                                                  function(x) unlist(strsplit(x, split="_"))[1])),
                      Feature = as.character(sapply(colnames(drugprof),
                                                    function(x) unlist(strsplit(x, split="_"))[2])))
rownames(colannot) = colnames(drugprof)


chan = c("#63b7af", "#035aa6", "#d8345f")
names(chan) = c("Cytoplasm", "Nuclei", "Cells")
fcol = c("#010059", "#107595", "#fdbfb3", 
         "#fcf594", "#aee7e8", "#ff9a76",
         "#8566aa")
names(fcol) = c("AreaShape", "Correlation", "Granularity",
                "Intensity", "Number", "RadialDistribution",
                "Texture")
ann_colors <- list(Level=chan,
                   Feature=fcol)

library(pheatmap)

paletteLength <- 200
heat_scale <- colorRampPalette(c("#78A1A7","white", "#C9788A"))(paletteLength)
maxval = 10
myBreaks <- c(seq(-maxval, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, maxval,
                  length.out=floor(paletteLength/2)))
ph = pheatmap(drugprof,
              # cluster_rows = row_clust$tree_row, 
              # cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)

save_pheatmap_pdf(ph, filename = paste0(figdir, 'Bray-big-heatmap.pdf'),
                  width = 10, height = 90)

drugsel = ph$tree_row$order[c(1:5, 6:11, 52:57, 263:265, 630:632, 728:740)]
subset_drugs = rownames(drugprof)[drugsel]


ph = pheatmap(drugprof[subset_drugs,],
              # cluster_rows = row_clust$tree_row, 
              # cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'Bray-selected-drugs.pdf'),
                  width = 10, height = 6)

drugprof[which(rowSums(abs(drugprof)) > quantile(rowSums(abs(drugprof)), 0.9)),]
top10 = drugprof[which(rowSums(abs(drugprof)) > quantile(rowSums(abs(drugprof)), 0.9)),]
ph = pheatmap(top10,
              # cluster_rows = row_clust$tree_row, 
              # cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_row = 0, 
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'Bray-top10percent-drugs.pdf'),
                  width = 10, height = 10)

ph = pheatmap(top10,
              # cluster_rows = row_clust$tree_row, 
              # cluster_cols = col_clust$tree_col,
              color = heat_scale,
              breaks = myBreaks, 
              show_colnames = F,
              annotation_col = colannot,
              annotation_colors = ann_colors,
              treeheight_col = 0)
save_pheatmap_pdf(ph, filename = paste0(figdir, 'Bray-newheatmap.pdf'),
                  width = 12, height = 10)

