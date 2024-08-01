setwd('G:/My Drive/lab files/Bryan Bergman/pancreatic islets analysis')
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(factoextra)
library(ggplot2)
library(enrichR)
library(FactoMineR)
library(reshape2)
library(biomaRt)

#PCA approach
setwd('G:/My Drive/lab files/Bryan Bergman/pancreatic islets analysis/')
cnts_mat = read.csv('melted_tpms unfiltered.csv')
head(cnts_mat)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attr = as.data.frame(listAttributes(ensembl))
ss1 = attr[grepl('symbol', attr$name),]
affyids <- as.vector(unique(cnts_mat$target_id))
conv_list = getBM(attributes = c('ensembl_transcript_id_version', 'hgnc_symbol'),
                  filters = 'ensembl_transcript_id_version',
                  values = affyids, 
                  mart = ensembl)
head(conv_list)

sample_table = read.csv('seq_sample_table.csv')
sample_table$treatment_time_hr
sample_table = sample_table[!sample_table$sample_ID==16,]
sample_table = sample_table[!sample_table$sample_ID==4,]
sample_table = sample_table[sample_table$treatment_time_hr==2,]
#or 16
cnts_mat$gene_symbol = conv_list$hgnc_symbol[match(cnts_mat$target_id, conv_list$ensembl_transcript_id_version)]
head(cnts_mat)

new_cnts = dcast(cnts_mat, sample_ID ~ gene_symbol, value.var = 'tpm', fun.aggregate = mean, na.rm=T)
row.names(new_cnts) = new_cnts$sample_ID
new_cnts$sample_ID=NULL

sample_table1 = sample_table
sample_table1$new_treatment_time = as.factor(sample_table1$treatment_time_hr)
sample_table1$new_sex = gsub('Male', 'M', sample_table1$Sex)
sample_table1$new_sex = gsub('Female', 'F', sample_table1$new_sex)
sample_table1$T2D_status=as.factor(sample_table1$T2D.)
sample_table1$sat_pat=as.factor(sample_table1$SAT.PAT)

sample_table1$SAT.PAT=NULL
sample_table1$T2D=NULL
sample_table1$treatment_time_hr=NULL
colnames(sample_table1) = gsub('new_', '', colnames(sample_table1))
sample_table1 = sample_table1[!is.na(sample_table1$sample_ID),]
row.names(sample_table1) = sample_table1$sample_ID
sample_table1$sample_ID = NULL
table(sample_table1$treatment_time)
new_cnts = new_cnts[row.names(new_cnts) %in%row.names(sample_table1), ]
summary(match(row.names(new_cnts), row.names(sample_table1)))

iris.pca <- PCA(new_cnts, graph = FALSE)

pca_groupings = function(group_variable) {
  new_sample_t = as.data.frame(sample_table[,colnames(sample_table) %in% group_variable])
  colnames(new_sample_t) = c('trait')
  new_sample_t$sample_ID = sample_table$sample_ID
  trait_cols =  unique(new_sample_t$trait)
  names(trait_cols) = MetBrewer::met.brewer('Moreau', length(trait_cols))
  new_sample_t$colors = names(trait_cols)[match(new_sample_t$trait, trait_cols)]
  indexing_cat = as.data.frame(row.names(new_cnts))
  colnames(indexing_cat) = c('sample_ID')
  indexing_cat$trait = new_sample_t$trait[match(indexing_cat$sample_ID, new_sample_t$sample_ID)]
  indexing_cat$pal = new_sample_t$colors[match(indexing_cat$sample_ID, new_sample_t$sample_ID)]
  
  g2 = fviz_pca_ind(iris.pca,
                    geom.ind = "point", # show points only (nbut not "text")
                    col.ind = indexing_cat$trait, # color by groups
                    palette = MetBrewer::met.brewer('Lakota', length(trait_cols)),
                    addEllipses = TRUE, # Concentration ellipses
                    legend.title = "Groups"
  ) + ggtitle(paste0('PCA plot grouping by ', gsub('new_', '', group_variable)))
  print(g2)
}
sample_table$new_BMI = ifelse(sample_table$BMI>30, 'above_30', 'below_30')
sample_table$new_treatment_time = as.factor(sample_table$treatment_time_hr)
sample_table$new_age = ifelse(sample_table$Age<50, '40_to_50', 'above_50')
sample_table$new_age = ifelse(sample_table$Age<40, 'below_40', paste0(sample_table$new_age))
sample_table$new_sex = gsub('Male', 'M', sample_table$Sex)
sample_table$new_sex = gsub('Female', 'F', sample_table$new_sex)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#21.00   38.75   47.50   45.75   53.25   62.00 
summary(sample_table$Age)

colnames(sample_table)
#[1] "sample_ID"         "HPAP.."            "SAT.PAT"           "Age"               "Sex"              
#[6] "BMI"               "Race.Ethnicity"    "T2D."              "treatment_time_hr"

g1 = pca_groupings('T2D.')
g2 = pca_groupings('new_BMI')
g3 = pca_groupings('new_age')
#g4 = pca_groupings('new_treatment_time')
g5 = pca_groupings('new_sex')
g6 = pca_groupings('SAT.PAT')
g7 = pca_groupings('Race.Ethnicity')
pdf(file = 'full PCA plots revised data.pdf')
gridExtra::grid.arrange(g1, g2, g5, g6, g7)
dev.off()
sample_table1 = sample_table
sample_table1$new_treatment_time = as.factor(sample_table1$treatment_time_hr)
sample_table1$new_sex = gsub('Male', 'M', sample_table1$Sex)
sample_table1$new_sex = gsub('Female', 'F', sample_table1$new_sex)
sample_table1$T2D_status=as.factor(sample_table1$T2D.)
sample_table1$sat_pat=as.factor(sample_table1$SAT.PAT)

sample_table1$SAT.PAT=NULL
sample_table1$T2D=NULL
sample_table1$treatment_time_hr=NULL
colnames(sample_table1) = gsub('new_', '', colnames(sample_table1))
sample_table1 = sample_table1[!is.na(sample_table1$sample_ID),]
row.names(sample_table1) = sample_table1$sample_ID
sample_table1$sample_ID = NULL
table(sample_table1$treatment_time)


######################################
#try adjustments
###set treatment time
time_pt = 2
  nn_sample = sample_table1[sample_table1$treatment_time==time_pt,]
  
  cnts_mat$log_cnts = log2(cnts_mat$est_counts + 1)
  new_cnts1 = dcast(cnts_mat, sample_ID ~ gene_symbol, value.var = 'log_cnts', fun.aggregate = mean, na.rm=T)
  new_cnts1 = new_cnts1[!is.na(new_cnts1$sample_ID),]
  
  row.names(new_cnts1) = new_cnts1$sample_ID
  
  new_cnts1 = new_cnts1[row.names(new_cnts1) %in% row.names(nn_sample),]
  new_cnts1$sample_ID=NULL
  new_cnts1[1:10,1:10]
  new_cnts1 = na.omit(new_cnts1)
  
  ####Filter by counts
  new_cnts1 = new_cnts1[,colSums(new_cnts1) >1]
  
  #see the grouping variable start with: new_sex, new_treatment_time, SAT.PAT, T2D.
  cnts_mat1 = new_cnts1
  cnts_mat1$dm = sample_table1$sat_pat[match(row.names(cnts_mat1), row.names(nn_sample))]
  library(limma)
  cnts_mat1 = cnts_mat1[!cnts_mat1$dm=='',]
  table(cnts_mat1$dm)
  design = model.matrix(~dm, data=cnts_mat1)
  head(design)
  table(cnts_mat1$dm)
  dim(design)
  remove_set = c('dm')
  new_cnts1 = as.data.frame(t(cnts_mat1[, !colnames(cnts_mat1) %in% remove_set]))
  fit = lmFit(new_cnts1, design)
  fit = eBayes(fit)
  row.names(fit)[1:10]
  
  res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH", sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
  head(res_table)
  
  new_cnts1[1:5,1:5]
  #change the name of the table
  test1 = as.data.frame(t(new_cnts1[row.names(new_cnts1)=='C19orf48P',]))
  test1$pS = sample_table1$sat_pat[match(row.names(test1), row.names(sample_table1))]
  
  write.csv(res_table, file = paste0(time_pt, 'hrs, results from limma on SAT over PAT - outlier removed.csv'), row.names = F)
  res1 = res_table
  res1$ID = row.names(res1)
  head(res1)
  
  #need to play around to assessing proper thresholds.  
  
  label_key = as.vector(res1$ID[res1$P.Value<0.01])
  res1$label2 = ifelse(res1$ID %in% label_key, paste0(res1$ID), '')
  table(res1$label2)
  res1$label_col1 = ifelse(res1$logFC>0, 'deepskyblue2', 'darkorange2')
  res1$label_col2 = ifelse(res1$P.Value<0.05, paste0(res1$label_col1), 'gray74')
  library(ggrepel)
  #Number of genes which will be labelled
  #Volcano plot
  #change labels here too
  pdf(file = paste0('Volcano Plot of ', time_pt, 'hrs PAT over SAT - age_sex adj.pdf'))
  ggplot(res1, aes(x=logFC, y=-log10(P.Value))) + theme_classic() +
    geom_point(aes(x=logFC, y=-log10(P.Value)), color=res1$label_col2) +
    geom_label_repel(aes(x=logFC, y=-log10(P.Value), label = res1$label2), color = res1$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, alpha = .6, segment.color = 'grey50')  +   ggtitle(paste0(time_pt, ' Volcano plot  PAT over SAT - age_sex adj'))
  
  dev.off()
  
  
  #####################tissue_enrichments
  library(clusterProfiler)
  library(enrichR)
  library(enrichplot)
  
  ff2 = res1[res1$P.Value < 0.3,]
  fc_dko = ff2$logFC
  names(fc_dko) <- ff2$ID
  fc_dko <- sort(fc_dko, decreasing = TRUE)
  fc_dko = fc_dko[!duplicated(names(fc_dko))]
  organism = "org.Hs.eg.db"
  ## Org.Hs.eg.db https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
  gse <-gseGO(
    geneList=fc_dko,
    ont = "ALL",
    OrgDb= organism,
    keyType = "SYMBOL",
    exponent = 1,
    minGSSize = 2,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH") 
  
  ## gsego https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/gseGO
  pdf(file = paste0(time_pt,' hrs GSEA output.pdf'))
  g2 = dotplot(gse, showCategory=10, split=".sign", color="pvalue") + facet_grid(.~.sign) + ggtitle('GSEA terms') + theme(axis.text.y=element_text(size=4))
  print(g2)
  dev.off()
  
  pdf(file = paste0(time_pt,  ' hrs GSEA Network graph.pdf'))
  x2<- pairwise_termsim(gse)
  g3 = emapplot(x2, showCategory = 10, color = "pvalue", cex_label_category=1)+ ggtitle('GSEA Network Graph') 
  print(g3)
  dev.off()
  
  pdf(file = paste0(time_pt,  ' hrs GSEA output and network.pdf'))
  g4 = gridExtra::grid.arrange(g2, g3, ncol=1)
  print(g4)
  dev.off()
  
  
  setEnrichrSite("Enrichr")
  dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")
  pp1 = res1[res1$P.Value < 0.05,]
  pp1 = pp1[pp1$logFC<1,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$ID
  
  #Change plot title names
  enriched <- enrichr(gg1, dbs1)
  g1 = plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('PAT upregulated genes - age_sex adj', names(enriched[1])))
  
  g2 = plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('PAT upregulated genes - age_sex adj', names(enriched[2])))
  
  
  g3 = plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('PAT upregulated genes - age_sex adj', names(enriched[3])))
  
  pdf(file = paste0(time_pt,'hrs PAT upregulated genes - age_sex adj.pdf'))
  gridExtra::grid.arrange(g1, g2, g3)
  dev.off()
  
  pp1 = res1[res1$P.Value < 0.05,]
  pp1 = pp1[pp1$logFC>1,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$ID
  enriched <- enrichr(gg1, dbs1)
  g1 = plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('PAT downregulated genes - age_sex adj', names(enriched[1])))
  
  g2 = plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('PAT downregulated genes - age_sex adj', names(enriched[2])))
  
  
  g3 = plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('PAT downregulated genes - age_sex adj', names(enriched[3])))
  
  pdf(file = paste0(time_pt, 'hrs PAT downregulated genes - age_sex adj.pdf'))
  gridExtra::grid.arrange(g1, g2, g3)
  dev.off()
  
  
  
  
  
run_full_de(2)


run_full_de(24)
run_full_de(2|24)





marekr_key = read.csv('cell marker set for heatmap.csv')


sample_table$heatID = paste0(sample_table$SAT.PAT, '-', sample_table$treatment_time_hr, 'hrs', sample_table$sample_ID )
markers_set = marekr_key$gene_symbol
cc2 = cnts_mat[cnts_mat$gene_symbol %in% markers_set,]
cc2$heatID = sample_table$heatID[match(cc2$sample_ID, sample_table$sample_ID)]
cc2$scaled_value  = ave(cc2$tpm, cc2$gene_symbol, FUN=scale)

cc3 = reshape2::dcast(cc2, heatID ~ gene_symbol, value.var = 'scaled_value', fun.aggregate = mean, na.rm=T)
cc3 = cc3[!is.na(cc3$heatID),]
row.names(cc3) = cc3$heatID
cc3$heatID=NULL

m2 = as.matrix(cc3)
m2[is.na(m2)] = 0
range(m2)
breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/1000)
m2 = na.omit(m2)
rownames(m2)


annotation_row = data.frame(
  rownames(m2)
)
rownames(annotation_row) = annotation_row$rownames.m2.
annotation_row$pat_sat = sample_table$SAT.PAT[match(rownames(annotation_row), sample_table$heatID)]
annotation_row$time = as.factor(sample_table$treatment_time_hr[match(rownames(annotation_row), sample_table$heatID)])


annotation_row$rownames.m2.=NULL
annotation_col = data.frame(
  colnames(m2)
)
rownames(annotation_col) = annotation_col$colnames.m2.
annotation_col$gene_category = marekr_key$category[match(annotation_col$colnames.m2., marekr_key$gene_symbol)]
annotation_col$colnames.m2.=NULL
pdf(file = 'Row zscore for select genes.pdf')
pheatmap(m2, color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(length(breaksList)), breaks = breaksList,  cluster_rows = T, annotation_row = annotation_row, annotation_col = annotation_col, cluster_cols = F, main='Row zscore for select traits', fontsize_row = 5, fontsize_col = 5)
dev.off()





