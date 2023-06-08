getwd()

require(dplyr)
require(stringr)

################################################################################
# IMPORT GO Database
################################################################################

Import_GO = function(min_size = 10, max_size = 500){
  require(dplyr)
  
  # Data Download
  if (!all(c('gene_association.mgi','go_terms.mgi') %in% list.files(getwd()))){
    system(command = 'wget https://www.informatics.jax.org/downloads/reports/gene_association.mgi.gz')
    system(command = 'gunzip gene_association.mgi.gz')
   
    system(command = 'wget https://www.informatics.jax.org/downloads/reports/go_terms.mgi')
  }
  
  GO = read.csv('gene_association.mgi', sep = '\t',skip = 3,header = F)
  GO = GO[c('V2','V5')]
  names(GO) = c('mgi_id','GO_term')
  
  GO_ID = read.csv('go_terms.mgi', sep = '\t', header = F)
  names(GO_ID) = c('Ont','GO_term','description')
  
  # Merging with description
  GO = left_join(GO, GO_ID, by = 'GO_term')
  
  # Filtering GO Term by Size
  GO_size = GO%>%group_by(GO_term)%>%summarise(Size=n())
  GO_term_filtered = GO_size[(GO_size$Size >= min_size & GO_size$Size <= max_size),][['GO_term']]
  GO_filtered = GO[GO$GO_term %in% GO_term_filtered,]

  return(GO_filtered)  
}

GO = Import_GO(min_size = 0, max_size = Inf)

################################################################################
# PREPROCESSING
################################################################################

# Data Loading 
cts = read.csv(file = 'read-counts.txt',sep = '\t',skip = 1)


# Filtering
cts = cts[cts[["CLIP.35L33G.bam"]] >= 30,]
cts = cts[cts[["RNA.control.bam"]] >= 30,]
cts = cts[cts[["RNA.siLin28a.bam"]] >= 30,]
cts = cts[cts[["RNA.siLuc.bam"]] >= 30,]
cts = cts[cts[["RPF.siLuc.bam"]] >= 80,]
cts = cts[cts[["RPF.siLin28a.bam"]] >= 80,]


# RPKM Normalization 
cts[["CLIP.35L33G.bam"]] = (as.numeric(cts[["CLIP.35L33G.bam"]]) /
                              ( as.numeric(cts[["Length"]])*sum(as.numeric(cts[["CLIP.35L33G.bam"]])) )) * 1000*1000000
cts[["RNA.control.bam"]] = (as.numeric(cts[["RNA.control.bam"]]) /
                              ( as.numeric(cts[["Length"]])*sum(as.numeric(cts[["RNA.control.bam"]])) )) * 1000*1000000
cts[["RNA.siLin28a.bam"]] = (as.numeric(cts[["RNA.siLin28a.bam"]]) /
                              ( as.numeric(cts[["Length"]])*sum(as.numeric(cts[["RNA.siLin28a.bam"]])) )) * 1000*1000000
cts[["RNA.siLuc.bam"]] = (as.numeric(cts[["RNA.siLuc.bam"]]) /
                              ( as.numeric(cts[["Length"]])*sum(as.numeric(cts[["RNA.siLuc.bam"]])) )) * 1000*1000000
cts[["RPF.siLuc.bam"]] = (as.numeric(cts[["RPF.siLuc.bam"]]) /
                              ( as.numeric(cts[["Length"]])*sum(as.numeric(cts[["RPF.siLuc.bam"]])) )) * 1000*1000000
cts[["RPF.siLin28a.bam"]] = (as.numeric(cts[["RPF.siLin28a.bam"]]) /
                            ( as.numeric(cts[["Length"]])*sum(as.numeric(cts[["RPF.siLin28a.bam"]])) )) * 1000*1000000


cts[['clip_enrichment']] = cts[['CLIP.35L33G.bam']]/ cts[['RNA.control.bam']]
cts[['rden_change']] = (cts[['RPF.siLin28a.bam']] / cts[['RNA.siLin28a.bam']]) / (cts[['RPF.siLuc.bam']] / cts[['RNA.siLuc.bam']])



################################################################################
# CONVERT GENE ID FOR GO ANALYSIS
################################################################################

# Converting Ensembl Gene ID to mgi ID
Convert_Ensembl_to_MGI = function(cts){
  
  # Loading biomaRt
  if (!('biomaRt.csv' %in% list.files(getwd()))){
    require(biomaRt)
    ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
    annot <- getBM(
      attributes = c(
        'mgi_id',
        'ensembl_gene_id_version'),
      filters = 'ensembl_gene_id_version',
      values = cts$Geneid,
      mart = ensembl)
    
    write.csv(annot, file = 'biomaRt.csv', row.names = F)
  } else {
    annot = read.csv('biomaRt.csv')
  }
  
  names(annot)[2] = 'Geneid'
  cts_merged = left_join(cts,annot,by='Geneid')
  
  # Sanity Check
  if (any(duplicated(cts_merged$Geneid))){
    message('Warning: Gene ID was duplicated when merging with biomaRt\n')
  }
  
  # Sanity Check
  if (any(is.na(cts_merged$mgi_id))){
    message(str_glue('Warning: A few Gene ID do not exit. {sum(is.na(cts_merged$mgi_id))} genes are removed. \n'))
    
    # Filter NA values
    cts_merged = cts_merged[!(is.na(cts_merged$mgi_id)),]
  }
  
  
  cts_merged$Geneid = cts_merged$mgi_id
  cts_merged$mgi_id = NULL
  
  return(cts_merged)
}


cts = Convert_Ensembl_to_MGI(cts)


################################################################################
# GO ANALYSIS: FUNCTION
################################################################################

# Filtering Fold Change Values From Count Matrix
FilterFC_FromCountMatrix <- function(gene_names, count_matrix, col_to_compare) {
  
  # Subsetting Count matrix based on Gene names
  rows_to_filter = (count_matrix$Geneid %in% gene_names)
  
  gene_set_FC <- count_matrix[rows_to_filter, ][[col_to_compare]]
  other_gene_FC <- count_matrix[!rows_to_filter, ][[col_to_compare]]
  
  # Return Results (list type)
  res = list(gene_set_FC, other_gene_FC)
  names(res) = c('gene_set_FC','other_gene_FC')
  return(res)
}

get_p_value = function(gene_set_counts, other_gene_counts){
  # Applying the Mann-Whitney U test
  p_value = wilcox.test(gene_set_counts, other_gene_counts, alternative = "two.sided")$p.value
  
  return(p_value)
}


get_lfc = function(gene_set_counts, other_gene_counts){
  # Get Log2 FoldChange
  FoldChange = (mean(gene_set_counts) / mean(other_gene_counts))
  LFC = log2(FoldChange)
  
  return(LFC)
}


################################################################################
# GO ANALYSIS: MAIN
################################################################################

list_FromGO = split(GO$mgi_id, GO$GO_term)

# Pregress Bar
require(progress)
pb <- progress_bar$new(
  format = " Progress: [:bar] :percent, Estimated time: :eta",
  total = length(names(list_FromGO))
)


results = list()
for (i in 1:length(names(list_FromGO))){
  pb$tick()

  GO_term = names(list_FromGO)[i]
  gene_names = list_FromGO[[GO_term]]

  # Sanity Check
  if (sum(cts$Geneid %in% gene_names) < 2){
    next
  }
  
  # Get Results
  results[['GO_term']][i] = GO_term
  results[['GS_size_all']][i] = length(gene_names)
  
  list_FC_rden = FilterFC_FromCountMatrix(gene_names, count_matrix = cts, col_to_compare = 'rden_change')
  list_FC_clip = FilterFC_FromCountMatrix(gene_names, count_matrix = cts, col_to_compare = 'clip_enrichment')
  
  results[['GS_size']][i] = length(list_FC_rden$gene_set_FC)
  results[['NotGS_size']][i] = length(list_FC_rden$other_gene_FC)
  
  results[['rden_p_value']][i] = get_p_value(list_FC_rden$gene_set_FC, list_FC_rden$other_gene_FC)
  results[['rden_enrichment']][i] = get_lfc(list_FC_rden$gene_set_FC, list_FC_rden$other_gene_FC)
  
  results[['clip_p_value']][i] = get_p_value(list_FC_clip$gene_set_FC, list_FC_clip$other_gene_FC)
  results[['clip_enrichment']][i] = get_lfc(list_FC_clip$gene_set_FC, list_FC_clip$other_gene_FC)
}

# Multiple Test Correction
results$rden_p_adj = p.adjust(results$rden_p_value, method = 'BH')
results$clip_p_adj = p.adjust(results$clip_p_value, method = 'BH')

# Remove NA
results = as.data.frame(results)
results = na.omit(results)

GO_subset = GO[c("GO_term", "Ont", "description")]
GO_subset = GO_subset[!duplicated(GO_subset),]

results_merged = left_join(results,GO_subset, by = 'GO_term')
write.csv(results_merged,'GO_results.csv', row.names = F)



################################################################################
# Visualization
################################################################################

# FUNCTION
visualize = function(results, filename = 'Figure_5_a.pdf'){

  require(ggplot2)
  require(ggrepel)
  require(colorspace)
  
  results_subset = results[results$clip_p_adj < 0.05 & results$rden_p_adj < 0.05,]
  results_subset$color = -log10(results_subset$clip_p_adj)
  results_subset[results_subset$color > 20,]$color= 20
  
  results_subset$label = paste0(results_subset$description," (",results_subset$GS_size, ") \n",
                                "C=",formatC(results_subset$clip_p_value, format = "e", digits = 1),
                                ", R=",formatC(results_subset$rden_p_value, format = "e", digits = 1))
  
  pathway_not_interest = c('mitochondrial inner membrane','nucleoplasm', 'cytosol','plasma membrane')
  #pathway_not_interest = c(pathway_not_interest,'plasma membrane','membrane')
  results_subset = results_subset[!(results_subset$description %in% pathway_not_interest),]
  
  pathway_of_interest = c('endoplasmic reticulum','endoplasmic reticulum lumen','endoplasmic reticulum membrane',
                          'Golgi apparatus', 'cell surface','mitochondrion','calcium ion binding',
                          'nucleus','cytoplasm','extracellular region','membrane')
  
  results
  
  p1 = ggplot(data = results_subset, aes(x = rden_enrichment, y = clip_enrichment, size = GS_size^2, color = color, label = label)) +
    geom_point(alpha = 0.7) +
    scale_size(range = c(0.7, 14)) +
    #  scale_color_gradient(low='yellow',high='red',limits=c(0, 20)) +
    #  scale_colour_gradientn(colours = c("darkred", "orange", "yellow"),breaks = c(0.1,0.01,0.001,0.0001), midpoint = 1.0e-10)+
    scale_color_continuous_sequential(palette = 'YlOrRd', limits=c(0, 20), breaks = seq(0,20,2)) +
    scale_x_continuous(breaks=seq(-3, 3, 0.5),limits = c(-3,3),expand = c(0,0)) +
    scale_y_continuous(breaks=seq(-1, 1.5, 0.5),limits = c(-1,1.5),expand = c(0,0)) +
    # geom_label_repel(size = 2, colour = 'black', 
    #                  data = results_subset%>%arrange(desc(GS_size),desc(color))%>%head, aes(label= label),
    #                  box.padding = 1) + 
    geom_label_repel(size = 2, colour = 'black',
                     data = results_subset[results_subset$description %in% pathway_of_interest,], aes(label= label),
                     box.padding = 0.8,
                     force = 10,
                     max.overlaps = Inf) +
    theme_bw() + 
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "black",size = 0.1,linetype = 'dashed'),
          panel.grid.minor = element_blank(),
          title         = element_text(size=11),
          axis.text.y   = element_text(size=8),
          axis.text.x   = element_text(size=8),
          axis.title.y  = element_text(size=9),
          axis.title.x  = element_text(size=10),
          legend.title  = element_text(size=10, angle = 90),
          legend.title.align = 0.5) +
    labs(title = 'Gene ontology term-enrichment analysis for CLIP and ribosome profiling',
         x = 'Enrichment level of LIN28A-bound CLIP tags (log2)',
         y = 'Ribosome density change upon Lin28a knockdown (log2)') +
    guides(size = 'none',
           color = guide_colourbar(title = "Enrichment Confidence: -log10(false discovery rate)",
                                   title.position = "left",
                                   barheight = 16)) +
    coord_cartesian(clip="off", xlim=c(-3, 3), ylim=c(-1, 1.5))
  p1
  ggsave(filename = filename,width = 20,height = 9,units = 'cm')
}


results = read.csv('GO_results.csv')
visualize(results)
