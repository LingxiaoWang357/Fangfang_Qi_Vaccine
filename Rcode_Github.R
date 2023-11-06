# Load packages -----------------------------------------------------------
Sys.setenv('R_MAX_VSIZE'=1024 * 100 * 1024^2)
set.seed(123)
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)


# Save R Image ------------------------------------------------------------
save.image("Rimage_pathway_231009.RData")


# load data ---------------------------------------------------------------
raw_mtx <- ReadMtx(
  mtx="../GSE171964_countsmatrix_v2.mtx.gz",
  cells = "../GSE171964_barcodes_v2.tsv.gz",
  features = "../GSE171964_feats_v2.tsv.gz",
  cell.column = 2,
  feature.column = 2,
  cell.sep = " ",
  feature.sep = " ",
  skip.cell = 1,
  skip.feature = 1,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)

Annotation_file <- as.data.frame(Annotation_file)
rownames(Annotation_file) <- Annotation_file[,1]

seurat_all <- CreateSeuratObject(counts = raw_mtx,
                                 meta.data = Annotation_file)

#Subset all PBMCs collected from Day 0,1,2,7
Idents(seurat_all) <- "day"
seurat_all_0127 <- subset(seurat_all, idents = c("0", "1", "2", "7"), downsample = 3000)
seurat_all_0127 <- NormalizeData(seurat_all_0127, scale.factor = 1e6)
seurat_all_0127 <- ScaleData(seurat_all_0127) 
seurat_all_0127 <- FindVariableFeatures(seurat_all_0127)
seurat_all_0127 <- RunPCA(seurat_all_0127)
seurat_all_0127 <- RunUMAP(seurat_all_0127, dims = 1:20)
Idents(seurat_all_0127) <- "clustnm"
DimPlot(seurat_all_0127, label = T, repel = T)

Idents(seurat_all_0127) <- factor(x = Idents(seurat_all_0127), 
                                  levels = c("C0_CD4 T",
                                             "C1_NK",
                                             "C2",
                                             "C3_CD14+ monocytes",
                                             "C4_CD16+ monocytes",
                                             "C5_B",
                                             "C6_CD8 T",
                                             "C7_cDC2",
                                             "C8_CD14+BDCA1+PD-L1+ cells",
                                             "C9_Platelets",
                                             "C10_Naive CD8 T",
                                             "C11_pDC",
                                             "C12_Tregs",
                                             "C13_cDC1",
                                             "C14_Plasmablasts",
                                             "C15_HPCs",
                                             "C16_NK T",
                                             "C17_Naive B"))
DimPlot(seurat_all_0127, label = T, repel = T, cols = c("#A47D4C",
                                                        "#818019",
                                                        "#414545",
                                                        "#253EF4",
                                                        "#006700",
                                                        "#CBCA8D",
                                                        "#F09737",
                                                        "#84308E",
                                                        "#FF999A",
                                                        "#3D9E68",
                                                        "#BF6757",
                                                        "#2B7FB7",
                                                        "#72FBFD",
                                                        "#936DC2",
                                                        "#87308A",
                                                        "#E9462A",
                                                        "#EB5AF9",
                                                        "#80F24B"))
DimPlot(seurat_all_0127,  cols = c("#A47D4C",
                                   "#818019",
                                   "#414545",
                                   "#253EF4",
                                   "#006700",
                                   "#CBCA8D",
                                   "#F09737",
                                   "#84308E",
                                   "#FF999A",
                                   "#3D9E68",
                                   "#BF6757",
                                   "#2B7FB7",
                                   "#72FBFD",
                                   "#936DC2",
                                   "#87308A",
                                   "#E9462A",
                                   "#EB5AF9",
                                   "#80F24B"))
DimPlot(seurat_all_0127, group.by = "day")
#Get gene expression Log(1+CPM)
All_cell_gene_count <- cbind(All_cell_gene_count,
                             as.data.frame(seurat_0127@meta.data$day),
                             as.data.frame(seurat_0127@meta.data$pt_id),
                             as.data.frame(seurat_0127@meta.data$clustnm))
colnames(All_cell_gene_count) <- c("cellid", "day", "patientid", "Cell_type")

gene_count <- as.data.frame(t(as.data.frame(seurat_0127@assays$RNA@data[c("IRF1", "STAT1", "JAK1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1"),])))

gene_count$cellid <- rownames(gene_count)

All_cell_gene_count <- left_join(All_cell_gene_count, gene_count, by = "cellid")
write.csv(All_cell_gene_count, "IFNG-related_logCPM.csv")


# functional study --------------------------------------------------------

# CD4, CD8, B cells, CD14+ monocytes, CD16+ monocytes, pDC, and cDC subsets.

IFN_features_large <- list(c("TANK","NR1H3","IFNL4","ADAR","IFITM3","IFITM2","CDC37","USP18","TREX1","TRIM6","OTOP1","DCST1","IFNLR1","TTLL12","YTHDF3","SAMHD1","LSM14A","IFNL2","IFNL3","IFNL1","TBK1","CNOT7","HCK","UBE2K","HPX","IFNE","STING1","IFI27","IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNAR1","IFNAR2","IFNB1","IRGM","IFNG","IFNGR1","IFNGR2","IFNW1","IL10RB","IRAK1","IRF1","IRF3","IRF7","JAK1","JAK2","ARG1","USP27X","MIR21","MMP12","MYD88","OAS1","OAS2","OAS3","YTHDF2","RBM47","PARP14","PPARG","MED1","METTL3","IFNK","MAVS","USP29","PTPN1","PTPN2","PTPN6","PTPN11","CACTIN","AZI2","SP100","STAT1","STAT2","TP53","TRAF3","TXK","TYK2","NR1H2","WNT5A","MUL1","ZBP1","TRIM56","PARP9","NLRC5","IFITM1","FADD","TRIM41","RNF185","ISG15","IKBKE","TBKBP1"))

# add module scores to all cells
seurat_all_0127 <- AddModuleScore(
  object = seurat_all_0127,
  features = IFN_features_large,
  name = 'IFN_Mediated_Pathways_large'
)

# export score of each cell to excel sheet
All_cell_module_score <- as.data.frame(seurat_all_0127@meta.data$barcode)

All_cell_module_score <- cbind(All_cell_module_score,
                               as.data.frame(seurat_all_0127@meta.data$day),
                               as.data.frame(seurat_all_0127@meta.data$pt_id),
                               as.data.frame(seurat_all_0127@meta.data$clustnm),
                               as.data.frame(seurat_all_0127@meta.data$IFN_Mediated_Pathways_large1))
colnames(All_cell_module_score) <- c("cellid", "day", "patientid", "Cell_type", "IFN_Mediated_Pathways_large")
head(All_cell_module_score)

write.csv(All_cell_module_score, "IFN_Mediated_Pathways_module_score_all_cells.csv")

# GSEA analysis of different cell types
unique(Idents(seurat_all_0127))
# CD4 T cell
seurat_all_0127_cd4tcell <- subset(seurat_all_0127, idents = "C0_CD4 T")
Idents(seurat_all_0127_cd4tcell) <- "day"
#DEGs
seurat_all_0127_cd4tcell_DEGs_1_0 <- FindMarkers(seurat_all_0127_cd4tcell, 
                                                 ident.1 = "1", 
                                                 ident.2 = "0",
                                                 min.pct = 0, 
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd4tcell_DEGs_1_0, "DEGs_cd4tcell_1_0.csv")

seurat_all_0127_cd4tcell_DEGs_2_0 <- FindMarkers(seurat_all_0127_cd4tcell, 
                                                 ident.1 = "2", 
                                                 ident.2 = "0",
                                                 min.pct = 0, 
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd4tcell_DEGs_2_0, "DEGs_cd4tcell_2_0.csv")

seurat_all_0127_cd4tcell_DEGs_7_0 <- FindMarkers(seurat_all_0127_cd4tcell, 
                                                 ident.1 = "7", 
                                                 ident.2 = "0",
                                                 min.pct = 0, 
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd4tcell_DEGs_7_0, "DEGs_cd4tcell_7_0.csv")


#GSEA
seurat_all_0127_cd4tcell_GSEA <- FindMarkers(seurat_all_0127_cd4tcell, 
                                             ident.1 = "1", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300,
                                             test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_cd4tcell_GSEA_FC <- seurat_all_0127_cd4tcell_GSEA$avg_log2FC
names(seurat_all_0127_cd4tcell_GSEA_FC) <- rownames(seurat_all_0127_cd4tcell_GSEA)
GSEA_seurat_all_0127_cd4tcell <- gseGO(geneList=seurat_all_0127_cd4tcell_GSEA_FC, 
                                       ont ="ALL", 
                                       keyType = "SYMBOL", 
                                       minGSSize = 5, 
                                       maxGSSize = 500, 
                                       pvalueCutoff = 0.05, 
                                       verbose = TRUE, 
                                       OrgDb = org.Hs.eg.db, 
                                       pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_cd4tcell, "GSEA_cd4tcell_day1vs0.csv")

#CD 8 T cell
seurat_all_0127_cd8tcell <- subset(seurat_all_0127, idents = "C6_CD8 T")
Idents(seurat_all_0127_cd8tcell) <- "day"

#DEGs
seurat_all_0127_cd8tcell_DEGs_1_0 <- FindMarkers(seurat_all_0127_cd8tcell, 
                                                 ident.1 = "1", 
                                                 ident.2 = "0",
                                                 min.pct = 0, 
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd8tcell_DEGs_1_0, "DEGs_cd8tcell_1_0.csv")

seurat_all_0127_cd8tcell_DEGs_2_0 <- FindMarkers(seurat_all_0127_cd8tcell, 
                                                 ident.1 = "2", 
                                                 ident.2 = "0",
                                                 min.pct = 0, 
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd8tcell_DEGs_2_0, "DEGs_cd8tcell_2_0.csv")

seurat_all_0127_cd8tcell_DEGs_7_0 <- FindMarkers(seurat_all_0127_cd8tcell, 
                                                 ident.1 = "7", 
                                                 ident.2 = "0",
                                                 min.pct = 0, 
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd8tcell_DEGs_7_0, "DEGs_cd8tcell_7_0.csv")

#GSEA
seurat_all_0127_cd8tcell_GSEA <- FindMarkers(seurat_all_0127_cd8tcell, 
                                             ident.1 = "1", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300,
                                             test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_cd8tcell_GSEA_FC <- seurat_all_0127_cd8tcell_GSEA$avg_log2FC
names(seurat_all_0127_cd8tcell_GSEA_FC) <- rownames(seurat_all_0127_cd8tcell_GSEA)
GSEA_seurat_all_0127_cd8tcell <- gseGO(geneList=seurat_all_0127_cd8tcell_GSEA_FC, 
                                       ont ="ALL", 
                                       keyType = "SYMBOL", 
                                       minGSSize = 5, 
                                       maxGSSize = 500, 
                                       pvalueCutoff = 0.05, 
                                       verbose = TRUE, 
                                       OrgDb = org.Hs.eg.db, 
                                       pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_cd8tcell, "GSEA_cd8tcell_day1vs0.csv")

# B cells
seurat_all_0127_bcell <- subset(seurat_all_0127, idents = "C5_B")
Idents(seurat_all_0127_bcell) <- "day"

#DEGs
seurat_all_0127_bcell_DEGs_1_0 <- FindMarkers(seurat_all_0127_bcell, 
                                              ident.1 = "1", 
                                              ident.2 = "0",
                                              min.pct = 0, 
                                              logfc.threshold = 0,
                                              only.pos = F,
                                              max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_bcell_DEGs_1_0, "DEGs_bcell_1_0.csv")

seurat_all_0127_bcell_DEGs_2_0 <- FindMarkers(seurat_all_0127_bcell, 
                                              ident.1 = "2", 
                                              ident.2 = "0",
                                              min.pct = 0, 
                                              logfc.threshold = 0,
                                              only.pos = F,
                                              max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_bcell_DEGs_2_0, "DEGs_bcell_2_0.csv")

seurat_all_0127_bcell_DEGs_7_0 <- FindMarkers(seurat_all_0127_bcell, 
                                              ident.1 = "7", 
                                              ident.2 = "0",
                                              min.pct = 0, 
                                              logfc.threshold = 0,
                                              only.pos = F,
                                              max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_bcell_DEGs_7_0, "DEGs_bcell_7_0.csv")

#GSEA

seurat_all_0127_bcell_GSEA <- FindMarkers(seurat_all_0127_bcell, 
                                          ident.1 = "7", 
                                          ident.2 = "0",
                                          min.pct = 0, 
                                          logfc.threshold = 0,
                                          only.pos = F,
                                          max.cells.per.ident = 300,
                                          test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_bcell_GSEA_FC <- seurat_all_0127_bcell_GSEA$avg_log2FC
names(seurat_all_0127_bcell_GSEA_FC) <- rownames(seurat_all_0127_bcell_GSEA)
GSEA_seurat_all_0127_bcell <- gseGO(geneList=seurat_all_0127_bcell_GSEA_FC, 
                                    ont ="ALL", 
                                    keyType = "SYMBOL", 
                                    minGSSize = 5, 
                                    maxGSSize = 500, 
                                    pvalueCutoff = 0.05, 
                                    verbose = TRUE, 
                                    OrgDb = org.Hs.eg.db, 
                                    pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_bcell, "GSEA_bcell_day7vs0.csv")

#cd14+ monocytes
seurat_all_0127_cd14monocytes <- subset(seurat_all_0127, idents = "C3_CD14+ monocytes")
Idents(seurat_all_0127_cd14monocytes) <- "day"

#DEGs
seurat_all_0127_cd14monocytes_DEGs_1_0 <- FindMarkers(seurat_all_0127_cd14monocytes, 
                                                      ident.1 = "1", 
                                                      ident.2 = "0",
                                                      min.pct = 0, 
                                                      logfc.threshold = 0,
                                                      only.pos = F,
                                                      max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd14monocytes_DEGs_1_0, "DEGs_cd14monocytes_1_0.csv")

seurat_all_0127_cd14monocytes_DEGs_2_0 <- FindMarkers(seurat_all_0127_cd14monocytes, 
                                                      ident.1 = "2", 
                                                      ident.2 = "0",
                                                      min.pct = 0, 
                                                      logfc.threshold = 0,
                                                      only.pos = F,
                                                      max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd14monocytes_DEGs_2_0, "DEGs_cd14monocytes_2_0.csv")

seurat_all_0127_cd14monocytes_DEGs_7_0 <- FindMarkers(seurat_all_0127_cd14monocytes, 
                                                      ident.1 = "7", 
                                                      ident.2 = "0",
                                                      min.pct = 0, 
                                                      logfc.threshold = 0,
                                                      only.pos = F,
                                                      max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd14monocytes_DEGs_7_0, "DEGs_cd14monocytes_7_0.csv")


#GSEA
seurat_all_0127_cd14monocytes_GSEA <- FindMarkers(seurat_all_0127_cd14monocytes, 
                                                  ident.1 = "2", 
                                                  ident.2 = "0",
                                                  min.pct = 0, 
                                                  logfc.threshold = 0,
                                                  only.pos = F,
                                                  max.cells.per.ident = 300,
                                                  test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_cd14monocytes_GSEA_FC <- seurat_all_0127_cd14monocytes_GSEA$avg_log2FC
names(seurat_all_0127_cd14monocytes_GSEA_FC) <- rownames(seurat_all_0127_cd14monocytes_GSEA)
GSEA_seurat_all_0127_cd14monocytes <- gseGO(geneList=seurat_all_0127_cd14monocytes_GSEA_FC, 
                                            ont ="ALL", 
                                            keyType = "SYMBOL", 
                                            minGSSize = 5, 
                                            maxGSSize = 500, 
                                            pvalueCutoff = 0.05, 
                                            verbose = TRUE, 
                                            OrgDb = org.Hs.eg.db, 
                                            pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_cd14monocytes, "GSEA_cd14monocytes_day2vs0.csv")


#cd16+ monocytes
seurat_all_0127_cd16monocytes <- subset(seurat_all_0127, idents = "C4_CD16+ monocytes")
Idents(seurat_all_0127_cd16monocytes) <- "day"

#DEGs
seurat_all_0127_cd16monocytes_DEGs_1_0 <- FindMarkers(seurat_all_0127_cd16monocytes, 
                                                      ident.1 = "1", 
                                                      ident.2 = "0",
                                                      min.pct = 0, 
                                                      logfc.threshold = 0,
                                                      only.pos = F,
                                                      max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd16monocytes_DEGs_1_0, "DEGs_cd16monocytes_1_0.csv")

seurat_all_0127_cd16monocytes_DEGs_2_0 <- FindMarkers(seurat_all_0127_cd16monocytes, 
                                                      ident.1 = "2", 
                                                      ident.2 = "0",
                                                      min.pct = 0, 
                                                      logfc.threshold = 0,
                                                      only.pos = F,
                                                      max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd16monocytes_DEGs_2_0, "DEGs_cd16monocytes_2_0.csv")

seurat_all_0127_cd16monocytes_DEGs_7_0 <- FindMarkers(seurat_all_0127_cd16monocytes, 
                                                      ident.1 = "7", 
                                                      ident.2 = "0",
                                                      min.pct = 0, 
                                                      logfc.threshold = 0,
                                                      only.pos = F,
                                                      max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cd16monocytes_DEGs_7_0, "DEGs_cd16monocytes_7_0.csv")

#GSEA

seurat_all_0127_cd16monocytes_GSEA <- FindMarkers(seurat_all_0127_cd16monocytes, 
                                                  ident.1 = "2", 
                                                  ident.2 = "0",
                                                  min.pct = 0, 
                                                  logfc.threshold = 0,
                                                  only.pos = F,
                                                  max.cells.per.ident = 300,
                                                  test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_cd16monocytes_GSEA_FC <- seurat_all_0127_cd16monocytes_GSEA$avg_log2FC
names(seurat_all_0127_cd16monocytes_GSEA_FC) <- rownames(seurat_all_0127_cd16monocytes_GSEA)
GSEA_seurat_all_0127_cd16monocytes <- gseGO(geneList=seurat_all_0127_cd16monocytes_GSEA_FC, 
                                            ont ="ALL", 
                                            keyType = "SYMBOL", 
                                            minGSSize = 5, 
                                            maxGSSize = 500, 
                                            pvalueCutoff = 0.05, 
                                            verbose = TRUE, 
                                            OrgDb = org.Hs.eg.db, 
                                            pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_cd16monocytes, "GSEA_cd16monocytes_day2vs0.csv")

#pDC
seurat_all_0127_pDc <- subset(seurat_all_0127, idents = "C11_pDC")
Idents(seurat_all_0127_pDc) <- "day"

#DEGs
seurat_all_0127_pDC_DEGs_1_0 <- FindMarkers(seurat_all_0127_pDc, 
                                            ident.1 = "1", 
                                            ident.2 = "0",
                                            min.pct = 0, 
                                            logfc.threshold = 0,
                                            only.pos = F,
                                            max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_pDC_DEGs_1_0, "DEGs_pDC_1_0.csv")

seurat_all_0127_pDC_DEGs_2_0 <- FindMarkers(seurat_all_0127_pDc, 
                                            ident.1 = "2", 
                                            ident.2 = "0",
                                            min.pct = 0, 
                                            logfc.threshold = 0,
                                            only.pos = F,
                                            max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_pDC_DEGs_2_0, "DEGs_pDC_2_0.csv")

seurat_all_0127_pDC_DEGs_7_0 <- FindMarkers(seurat_all_0127_pDc, 
                                            ident.1 = "7", 
                                            ident.2 = "0",
                                            min.pct = 0, 
                                            logfc.threshold = 0,
                                            only.pos = F,
                                            max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_pDC_DEGs_7_0, "DEGs_pDC_7_0.csv")

#GSEA

seurat_all_0127_pDc_GSEA <- FindMarkers(seurat_all_0127_pDc, 
                                        ident.1 = "2", 
                                        ident.2 = "0",
                                        min.pct = 0, 
                                        logfc.threshold = 0,
                                        only.pos = F,
                                        max.cells.per.ident = 300,
                                        test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_pDc_GSEA_FC <- seurat_all_0127_pDc_GSEA$avg_log2FC
names(seurat_all_0127_pDc_GSEA_FC) <- rownames(seurat_all_0127_pDc_GSEA)
GSEA_seurat_all_0127_pDc <- gseGO(geneList=seurat_all_0127_pDc_GSEA_FC, 
                                  ont ="ALL", 
                                  keyType = "SYMBOL", 
                                  minGSSize = 5, 
                                  maxGSSize = 500, 
                                  pvalueCutoff = 0.05, 
                                  verbose = TRUE, 
                                  OrgDb = org.Hs.eg.db, 
                                  pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_pDc, "GSEA_pDC_day2vs0.csv")

#cDC1
seurat_all_0127_cDC1 <- subset(seurat_all_0127, idents = "C13_cDC1")
Idents(seurat_all_0127_cDC1) <- "day"

#DEGs
seurat_all_0127_cDC1_DEGs_1_0 <- FindMarkers(seurat_all_0127_cDC1, 
                                             ident.1 = "1", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cDC1_DEGs_1_0, "DEGs_cDC1_1_0.csv")

seurat_all_0127_cDC1_DEGs_2_0 <- FindMarkers(seurat_all_0127_cDC1, 
                                             ident.1 = "2", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cDC1_DEGs_2_0, "DEGs_cDC1_2_0.csv")

seurat_all_0127_cDC1_DEGs_7_0 <- FindMarkers(seurat_all_0127_cDC1, 
                                             ident.1 = "7", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cDC1_DEGs_7_0, "DEGs_cDC1_7_0.csv")


#GSEA

seurat_all_0127_cDC1_GSEA <- FindMarkers(seurat_all_0127_cDC1, 
                                         ident.1 = "2", 
                                         ident.2 = "0",
                                         min.pct = 0, 
                                         logfc.threshold = 0,
                                         only.pos = F,
                                         max.cells.per.ident = 300,
                                         test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_cDC1_GSEA_FC <- seurat_all_0127_cDC1_GSEA$avg_log2FC
names(seurat_all_0127_cDC1_GSEA_FC) <- rownames(seurat_all_0127_cDC1_GSEA)
GSEA_seurat_all_0127_cDC1 <- gseGO(geneList=seurat_all_0127_cDC1_GSEA_FC, 
                                   ont ="ALL", 
                                   keyType = "SYMBOL", 
                                   minGSSize = 5, 
                                   maxGSSize = 500, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE, 
                                   OrgDb = org.Hs.eg.db, 
                                   pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_cDC1, "GSEA_cDC1_day2vs0.csv")

#cDC2
seurat_all_0127_cDC2 <- subset(seurat_all_0127, idents = "C7_cDC2")
Idents(seurat_all_0127_cDC2) <- "day"

#DEGs
seurat_all_0127_cDC2_DEGs_1_0 <- FindMarkers(seurat_all_0127_cDC2, 
                                             ident.1 = "1", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cDC2_DEGs_1_0, "DEGs_cDC2_1_0.csv")

seurat_all_0127_cDC2_DEGs_2_0 <- FindMarkers(seurat_all_0127_cDC2, 
                                             ident.1 = "2", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cDC2_DEGs_2_0, "DEGs_cDC2_2_0.csv")

seurat_all_0127_cDC2_DEGs_7_0 <- FindMarkers(seurat_all_0127_cDC2, 
                                             ident.1 = "7", 
                                             ident.2 = "0",
                                             min.pct = 0, 
                                             logfc.threshold = 0,
                                             only.pos = F,
                                             max.cells.per.ident = 300) %>% arrange(desc(avg_log2FC))
write.csv(seurat_all_0127_cDC2_DEGs_7_0, "DEGs_cDC2_7_0.csv")

#GSEA
seurat_all_0127_cDC2_GSEA <- FindMarkers(seurat_all_0127_cDC2, 
                                         ident.1 = "2", 
                                         ident.2 = "0",
                                         min.pct = 0, 
                                         logfc.threshold = 0,
                                         only.pos = F,
                                         max.cells.per.ident = 300,
                                         test.use = "DESeq2") %>% arrange(desc(avg_log2FC))
seurat_all_0127_cDC2_GSEA_FC <- seurat_all_0127_cDC2_GSEA$avg_log2FC
names(seurat_all_0127_cDC2_GSEA_FC) <- rownames(seurat_all_0127_cDC2_GSEA)
GSEA_seurat_all_0127_cDC2 <- gseGO(geneList=seurat_all_0127_cDC2_GSEA_FC, 
                                   ont ="ALL", 
                                   keyType = "SYMBOL", 
                                   minGSSize = 5, 
                                   maxGSSize = 500, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE, 
                                   OrgDb = org.Hs.eg.db, 
                                   pAdjustMethod = "BH")
write.csv(GSEA_seurat_all_0127_cDC2, "GSEA_cDC2_day2vs0.csv")





