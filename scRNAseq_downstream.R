### Anna Mathioudaki ### 
####### 2020 Feb #######

# Downstream analysis of 2 scRNA-seq datasets after integration
# Using Seurat and Monocle 
# Cells: Mesenchymal Stromal Cells
# Different individuals were mixed in one 10x Genomics run
# Prior donor demultiplexing using cellSNP and vireo

# Load required packages
pacman::p_load(devtools, Seurat, Matrix, tidyverse, zinbwave,
               clusterProfiler, org.Hs.eg.db, umap,
               pals, plotly, scater, cowplot, ClusterMap, hash, 
               monocle, clues, sctransform, pheatmap, RColorBrewer)

# Create parser argument
parser <- ArgumentParser()
parser$add_argument("-tenxrun1", "--tenxrun1", required=TRUE, help="Full path for Cell Ranger output folder -filtered_feature_bc_matrix- for first run")
parser$add_argument("-projectname1", "--projectname1", required=TRUE, help="Project name for seurat object and cell ids")
parser$add_argument("-vireo1", "vireo1", required=TRUE, help="donor_ids.tsv output of vireo that indicates demultiplexing of mixed individuals in single cell experiment")

parser$add_argument("-10xrun2", "--10xrun2", required=TRUE, help="Full path for Cell Ranger output folder -filtered_feature_bc_matrix- for second run")
parser$add_argument("-projectname2", "--projectname2", required=TRUE, help="Project name for seurat object and cell ids")
parser$add_argument("-vireo2", "vireo2", required=TRUE, help="donor_ids.tsv output of vireo that indicates demultiplexing of mixed individuals in single cell experiment")

parser$add_argument("-projectname", "--projectname", required=TRUE, help="Name for the project when loaded to seurat")
parser$add_argument("-n_pcs", "-n_pcs", required=TRUE, help="Number of principal components (pcs) to be used for further downstream analysis as determined by elbow plot.")
parser$add_argument("-mindist", "--mindist", required=TRUE, help="min_dist argument used for umap and clustering - after optimization")
parser$add_argument("-nclusters", "--nclusters", required=TRUE, help="Number of clusters identified from clustering")
parser$add_argument("-dims", "--dims", required=TRUE, help="Number of dimensions to include in monocle clustering")

args <- parser$parse_args()

######################################################
# Define functions 
filter_mitochondrial_genes <- function(seurat_object){
  mt_genes <- rownames(seurat_object)[grep("^MT-",rownames(seurat_object))]
  length(mt_genes)
  C <- GetAssayData(object = seurat_object, slot = "counts") 
  if (length(mt_genes) > 0 ){
    cat(length(mt_genes), "genes were identified as mitochondrial. This information will be added in the metadata.")
    percent_mito <- colSums(C[mt_genes, ])/Matrix::colSums(C)*100
    seurat_object <- AddMetaData(seurat_object, percent_mito, col.name = "percent_mito")
    return(seurat_object)
  }
  else{
    print("Zero mitochondrial genes were detected", quote=FALSE)
  }
} 

combined_seurat_final <- function(seurat1_matrix,
                                  seurat2_matrix,
                                  min_cells, 
                                  nfeat_min,
                                  seurat1name, 
                                  seurat2name,
                                  project){
  # Create Seurat object
  sample1_seurat <- CreateSeuratObject(counts=seurat1_matrix, min.cells=min_cells, project=seurat1name) # min.cells features presented in at least min.cells cells
  sample2_seurat <- CreateSeuratObject(counts=seurat2_matrix, min.cells=min_cells, project=seurat2name) # min.cells features presented in at least min.cells cells
  combined_seurat <- merge(sample1_seurat, y = sample2_seurat, add.cell.ids = c(seurat1name, seurat2name), project = project)
  
  # store mitochondrial percentage in object meta data
  combined_seurat <- PercentageFeatureSet(combined_seurat, pattern = "^MT-", col.name = "percent.mt")
  # Plot features for seurat object. This is a QC representation 
  p <- VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
  q <- FeatureScatter(object = combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")   # check for the scatterplot of N features vs N coutns RNA
  # Filter seurat object for cells with low gene detection (cutoffs determined by the chemistry)
  #nFeature_RNA -->  # of detected genes
  #did not filter for high gene detection and thus douplets, will do based on genotype (vireo output)
  a <<- as.numeric(nfeat_min)
  combined_seurat <- subset(x = combined_seurat, subset = nFeature_RNA > a & percent.mt < 10)
  combined_seurat <- SCTransform(combined_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
  r <- VariableFeaturePlot(object = combined_seurat)
  
  final_seurat <<- combined_seurat
  c(plot(p), plot(q), plot(r))
  
  VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
}

get_marker_genes <- function(seurat_object, number_of_clusters){
  #seurat_object: definition of the object marker-genes per cluster should be idintified
  #number_of_clusters: number of clusters occuring in the particular object after clustering algorithms
  genes_per_cluster <- hash() # this is where u store the genes per cluster 
  for (i in 0:as.numeric(number_of_clusters)){
    cluster_markers <- FindMarkers(object=seurat_object, ident.1=i, min.pct = 0.25)
    #min.pct argument requires a feature to be detected at minimum percentage in either of the 2 groups of cells
    genes_per_cluster[[paster(c("cluster_", i), collapse="")]] <- as.list(cluster_genes)
  }
  print(genes_per_cluster)
}

seurat_to_monocle <- function(seurat_object){
  data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  #Construct monocle cds
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  cds <- monocle_cds
  # Estimate size factors and dispersions
  cds <<- estimateSizeFactors(cds)
  cds <<- estimateDispersions(cds)
  #http://cole-trapnell-lab.github.io/monocle-release/docs/
  #FPKM/TPM values are generally log-normally distributed
  #UMIs or read counts are better modeled with the negative binomial / do not normalize prior to creating CellDataSet
}

######################################################
#Create seurat object from the 10x matrices
sample1matrix <- Read10X(data.dir=args$tenxrun1, gene.column = 2, unique.features = TRUE)
sample2matrix <- Read10X(data.dir=args$tenxrun1, gene.column = 2, unique.features = TRUE)
combined_seurat_final(sample1matrix, sample2matrix, 2, 500, args$projectname1, args$projectname2, args$projectname)

######################################################
#Dimensionality reduction and clustering
# run PCA 
final_seurat <- RunPCA(object = final_seurat, verbose = FALSE)#, pc.genes = final_seurat@var.genes)

pdf(paste0(args$projectname, "_PCA.pdf"), width=10, height=10)
DimPlot(object = final_seurat)
dev.off()
# stastically significant PCs 
pdf(paste0(args$projectname, "_elbow.pdf"), width=10, height=10)
ElbowPlot(object = final_seurat)
dev.ofF()

#Run tSNE
final_seurat  <- RunTSNE(object = final_seurat, dims = 1:30)

#Run UMAP
#First run for different min.dist values to see how the clusters seperate
for (min_dist in c(0.001, 0.01, 0.05, 0.1, 0.5, 1)){
  testSeurat_0 <- RunUMAP(object = final_seurat, reduction = "pca", dims = 1:30, min.dist = min_dist)
  testSeurat_0 <- FindNeighbors(object = testSeurat_0, dims=1:30, verbose=FALSE)
  testSeurat_0 <- FindClusters(object = testSeurat_0, algorithm=1, n.start=100, 
                               n.iter=200, pc.use = 1:args$n_pcs, resolution=0.5, 
                               save.SNN=T, do.sparse=T, verbose=FALSE)
  pdf(paste0("testumap_mindist", min_dist, ".pdf"))
  p1 <- DimPlot(testSeurat_0, reduction = "umap") + 
    ggtitle(paste0("min_dist=", min_dist, ", 1st PC used "))
  plot(p1)
  dev.off()
}

# Run UMAP
final_seurat <- RunUMAP(object = final_seurat, reduction = "pca", dims = 1:30, min.dist = args$mindist)

# Clustering
final_seurat <- FindNeighbors(object = final_seurat, dims=1:30, verbose=FALSE)
final_seurat <- FindClusters(object = final_seurat,  algorithm = 1,
                              n.start = 100, n.iter = 200 ,pc.use = 1:args$n_pcs, 
                              resolution = 0.5, save.SNN = T, do.sparse = T)

# Run & plot tSNE
pdf(paste0(args$projectname, "_tsne.pdf"), width=10, height=10)
DimPlot(final_seurat, reduction = "tsne")
dev.off()

pdf(paste0(args$projectname, "_umap.pdf"), width=10, height=10)
DimPlot(final_seurat, reduction = "umap")
dev.off()

######################################################
#Get an overview of which dataset each cell comes from
pdf("umap_dataset.pdf")
DimPlot(final_seurat, group.by = "orig.ident", reduction="umap", combine = F)
dev.off()
cat("Origin of each cell per cluster\n", table(final_seurat@active.ident, final_seurat@meta.data$orig.ident))


######################################################
# Now lets see in how many individuals each of the cluster corresponds
# Vireo
sample_vireo <- data.frame()
for (file in c(args$vireo1, args$vireo2)){
  vireo <- read.table(file, header=T)
  sample_vireo_i <- data.frame(vireo$cell,vireo$best_singlet)
  colnames(sample_vireo_i) <- c('cell', 'donor')
  tmp_cell <- c()
  for (i in sample_vireo_i$cell){
    cell <- strsplit(i, split="-")[[1]][1]
    tmp_cell <- c(tmp_cell, cell)
  }
  sample_vireo_i$cell <- tmp_cell
  sample_vireo <- rbind(sample_vireo, sample_vireo_i)
  rm(tmp_cell)
  rm(sample_vireo_i)
}

tmp <- data.frame()
tmp_dataset <- c()
tmp_cell <- c()
tmp_cluster <- c()

cells_not_in_vireo <- c()
for (i in names(final_seurat@active.ident)){
  j <- strsplit(i, split="_")[[1]][2] #cell
  if (j %in% sample_vireo$cell){
    tmp_dataset <- c(tmp_dataset, strsplit(i, split="_")[[1]][1])
    tmp_cell <- c(tmp_cell, j)
    tmp_cluster <- c(tmp_cluster, final_seurat@active.ident[[i]])  
  }
  else {
    cells_not_in_vireo <- c(cells_not_in_vireo, j) #get cells that are not in the vireo dataset
  }
}

tmp <- cbind(tmp_dataset, tmp_cell, tmp_cluster)
rm(tmp_dataset, tmp_cell, tmp_cluster)
colnames(tmp) <- c('dataset', 'cell', 'cluster')
tmp <- data.frame(tmp)

tmp_donor <- c()
for (i in tmp$cell){
  d <- sample_vireo[which(sample_vireo$cell==i), ]$donor
  d <- as.character(d)[1]
  tmp_donor <- c(tmp_donor, d)
}

final <- cbind(tmp, tmp_donor)
colnames(final) <- c('dataset', 'cell', 'cluster', 'donor')

print(table(final$donor, final$cluster, final$dataset))


######################################################
# Functional profiles per cluster using clusterProfiler
# MF -> Molecular function
# CC -> Cellular compound
# BP -> Biological processes

#Background dataset --> all genes expressed in each cluster
#GO for each cluster using ClusterProfiler
background1 <- AverageExpression(object=final_seurat)
colnames(background1$RNA) <- paste0('X', colnames(background1$RNA))

for (i in 0:args$n_clusters){#number of clusters / clusters start from 0 
  tmp <- data.frame(data.frame(background1$RNA)[ ,paste0('X', as.numeric(i))], row.names = row.names(background1$RNA))
  colnames(tmp) <- 'cluster'
  #background <- as.list(rownames(background1))
  background <- rownames(tmp[which(tmp$cluster!=0), , drop = FALSE]) #this is background based on all the genes expressed in this cluster
  markers <- FindMarkers(object = final_seurat, ident.1 = i, min.pct = 0.25) ######### Convert to ensembl ids
  for (go in c("MF", "CC", "BP")){
    go_enrich <- enrichGO(gene = rownames(markers),
                          OrgDb = 'org.Hs.eg.db',
                          keyType = "SYMBOL",
                          ont = go,
                          universe = as.character(background))
    pdf(paste0(args$projectname, "cluster", i, "_", go, ".pdf"), width=11, height=10)
    plot(dotplot(go_enrich, showCategory = 20, title = c(paste0("cluster_", i, " ", go), paste(c(" ", go), collapse = "")), font.size=12))
    dev.off()
  }
}

######################################################
# Identify marker genes per cluster
#Cells with a value > 0 represent cells with expression above the population mean 
#(a value of 1 would represent cells with expression 1SD away from the population mean).

msc_markers <- FindAllMarkers(final_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # logFC maybe not needed and might lose information 
msc_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(msc_markers, file=paste0(args$projectmame, "_markers.tsv"), quote=FALSE, sep="\t", row.names = T)

top10 <- msc_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) # look at everything

pdf(paste0(args$projectname, "_topmarkergenes.pdf"))
DoHeatmap(final_seurat, features = top10$gene, draw.line=TRUE) + theme(axis.text.y=element_text(size = 6))
dev.off()

######################################################
# Load seurat object to monocle for trajectory analysis 
#Extract data, phenotype data, and feature data from the SeuratObject
seurat_to_monocle(final_seurat)
cds <- detectGenes(cds, min_expr = 0.1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

#Count number of cells genes important for the MSC lineage are expressed
int_genes <- c("THY1", "CD44", "ENG", "VCAM1", "ALCAM", "NT5E", "PPARG", "PTPRC", "CD34", "HLA-DRA", 
               "CD19", "CD14", "RUNX2", "TGFBR2", "SMAD1", "GSK3B", "TOB1", "LRP6", "NEAT1", "NR2F2", 
               "HDAC9", "ALPL", "COL1A1", "FABP4", "GAPDH", "NES", "ADIPOQ", "LEP")
numcells_pergene <- fData(cds)[int_genes, ]
numcells_pergene <- numcells_pergene[order(numcells_pergene$num_cells_expressed), ]

print(numcells_pergene)

######################################################
#r Monocle-Dimensionality reduction
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

# setOrderingFilter function marks genes used for clustering in subsequent calls to clusterCells
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
# plot_ordering_genes function shows how variability (dispersion) in a gene's expression depends on the average expression across cells
#red line shows Monocle's expectation of the dispersion based on this relationship. 
#The genes we marked for use in clustering are shown as black dots, while the others are shown as grey dots.
pdf("dispersion_monocle.pdf")
plot_ordering_genes(cds)
dev.off()

#table(disp_table$mean_expression>=0.1)

# PCA dimensionality reduction for achieving shorter running times
pdf("pca_variance_monocle.pdf")
plot_pc_variance_explained(cds, return_all = FALSE)
dev.off()

#In order to compare the clustering results, we are usign Adjusted Rand Index (ARI) from the clues package
#Well explained: https://davetang.org/muse/2017/09/28/rand-index-versus-adjusted-rand-index/
#The ARI (HA) suggests that the clustering of cells is quite different when including different numbers of principal components. 

cluster_dim_eval <- list()
dimensions_tested <- list()

for (i in 3:10){
  dim <- paste0("dim_", i) 
  dimensions_tested <- c(dimensions_tested, dim)
  cds <- reduceDimension(cds, max_components = 2, num_dim = as.numeric(i),
                            reduction_method = 'tSNE', verbose = T)
  cds <- clusterCells(cds, num_clusters = 15) # requesting 1 to 10 clusters
  plot(plot_cell_clusters(cds, 1, 2) + labs(main=dim))
  cluster_dim_eval[[dim]] <- pData(cds)$Cluster
}


# In order to compare the clustering results, we are usign Adjusted Rand Index (ARI) from the clues package
pairwise_ARI <- data.frame()
for (i in dimensions_tested){
  for (j in dimensions_tested){
    pairwise_ARI[i, j] <- as.numeric(adjustedRand(as.numeric(cluster_dim_eval[[i]]), as.numeric(cluster_dim_eval[[j]]), randMethod="HA"))
  }
}    
pdf("ARI.pdf")
pheatmap(as.matrix(pairwise_ARI), col=brewer.brbg(9), cluster_cols = F, cluster_rows = F, main="Adjusted random index-PCs for clustering")
dev.off()

#Include dimensions based on the ARI output
cds <- reduceDimension(cds, max_components = 2, num_dim = args$dims,
                          reduction_method = 'tSNE', verbose = TRUE)
cds <- clusterCells(cds, num_clusters = 15)

######################################################
# Trajectory analysis 
# Use on of the colnames(pData(cds)) inside the fullModelFormulaStr
diff_timepoint <- differentialGeneTest(cds, relative_expr = FALSE, fullModelFormulaStr="~Cluster") #Runtime ~20min
ordering_sites <- row.names(subset(diff_timepoint, qval < 0.01))[1:1000] # qvalue cutoff -> 0.01 to 0.1
cds <- setOrderingFilter(cds, ordering_sites)

cds <- reduceDimension(cds, max_components = 2,
                           residualModelFormulaStr="~as.numeric(num_genes_expressed)",
                           reduction_method = 'DDRTree')
cds <- orderCells(cds)

pdf("trajectory_per_cluster.pdf")
plot_cell_trajectory(cds, color_by = "State")
dev.off()

# Trajectory for top10 genes per cluster
for (i in 0:args$n_clusters){
  range <- (1+(i*10)):((i*10)+10)
  print(top10$gene[range])
  pdf(paste0("trajectory_topgenescluster", i, ".pdf"))
  plot(plot_cell_trajectory(cds, markers = c(top10$gene[range]), use_color_gradient = TRUE)) + ggplot2::ggtitle(paste0("Trajectory plot for top10 genes for cluster", i+1))
  dev.off()
}

#plot_cell_trajectory(cds, markers = c("THY1", "VCAM1", "NEAT1", "NR2F2", "LRP6", "GSK3B", "NES", "ALPL", "RUNX2"), use_color_gradient = TRUE)
