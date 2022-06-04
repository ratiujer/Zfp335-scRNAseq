library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(Nebulosa)
library(EnhancedVolcano)
#read in count matrix for each of the conditions of interest and create Seurat Objs
setwd("G:/DN4-10X/counts")

## Names of folders containing the expression matrices for each sample
con<-Read10X("Control_raw_feature_bc_matrix")
zfp<-Read10X("Zfp335-cKO_raw_feature_bc_matrix")

con<-CreateSeuratObject(con,project="Con-dn4",min.cells=3, min.features = 500)
zfp<-CreateSeuratObject(zfp,project="Zfp335_cKO-dn4",min.cells=3, min.features = 500)
saveRDS(con,'con_seurat.rds')
saveRDS(zfp,'zfp_seurat.rds')

###Reprocess data - load con and zfp seurat objects

zfp<-readRDS("zfp_seurat.rds")
con<-readRDS("con_seurat.rds")
MitochondrialGENES <- readRDS("G:/DN4-10X/counts/MitochondrialGENES.rds")
mouse_cell_cycle_genes <- readRDS("G:/DN4-10X/counts/mouse_cell_cycle_genes.rds")


MitoPctC <- Matrix::colSums(con@assays$RNA@counts[MitochondrialGENES, ])/Matrix::colSums(con@assays$RNA@counts)
con <- AddMetaData(object = con, metadata = MitoPctC , col.name = "mito")

MitoPctZ <- Matrix::colSums(zfp@assays$RNA@counts[MitochondrialGENES, ])/Matrix::colSums(zfp@assays$RNA@counts)
zfp <- AddMetaData(object = zfp, metadata = MitoPctZ , col.name = "mito")

con[["geno"]]<-'WT'
zfp[["geno"]]<-'Zfp335cKO'

con<-subset(con,subset = nFeature_RNA > 1000 & nCount_RNA > 1000 & mito < 0.075)
zfp<-subset(zfp,subset = nFeature_RNA > 1000 & nCount_RNA > 1000 & mito < 0.075)

reference.list <- c(con, zfp)
reference.list<-lapply(X=reference.list, FUN=function(x) {
  x<- NormalizeData(x)
  x<-FindVariableFeatures(x,selection.method="vst", nfeatures=6000)
})

dn4.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 6000)

dn4.integrated<-IntegrateData(anchorset=dn4.anchors,dims=1:30)


dn4.integrated <- ScaleData(dn4.integrated,verbose=FALSE)

g2m.genes<-mouse_cell_cycle_genes$g2m.genes
s.genes<-mouse_cell_cycle_genes$s.genes

dn4.integrated<-CellCycleScoring(dn4.integrated,s.features = s.genes, g2m.features = g2m.genes)

dn4.int.reg<-ScaleData(dn4.integrated,vars.to.regress = c("S.Score","G2M.Score"),features=rownames(dn4.integrated))
dn4.int.reg<-RunPCA(dn4.int.reg,features=VariableFeatures(dn4.int.reg),nfeatures.print = 10)
dn4.int.reg<-JackStraw(dn4.int.reg,num.replicate = 100,dims = 40)
dn4.int.reg<-ScoreJackStraw(dn4.int.reg,dims=1:40)

#dn4.int.reg[['orig.cluster']]<-Idents(dn4.int.reg)
dn4.int.reg<-FindNeighbors(dn4.int.reg,dims=1:40)
dn4.int.reg<-FindClusters(dn4.int.reg,resolution = 0.5)
dn4.int.reg<-RunUMAP(dn4.int.reg,dims=1:40)
DimPlot(dn4.int.reg, label=T, label.size = 5, label.box = T)

new.ids<-c('DN4_1','DN4_2','Mat_1','DN4_3','Mat_2','DN4_4','gd17','Mat_3','gd1','DN4_5')
names(new.ids)<-levels(dn4.int.reg)
dn4.int.reg<-RenameIdents(dn4.int.reg,new.ids)
dn4.int.reg[['cellType']]<-Idents(dn4.int.reg)

markers<-FindAllMarkers(dn4.int.reg,assay = 'RNA')

DotPlot(dn4.int.reg,assay = 'RNA',c('Ptcra','Mki67','Cdk1','Rag1','Cd8b1','Cd8a','Cd4','Trbc1','Trac','Trdc','Tcrg-C1','Tcrg-C2','Tcrg-C4','Sox13','Rorc','Maf','Nkg7','Il2rb','S1pr1','Il7r'), cols = 'RdYlBu', col.max = 0.8,scale.by = 'radius',cluster.idents = T, dot.scale = 10, dot.min = 0.1)+RotatedAxis()

## Remove maturing ab and gd T cells to analyze only DN4 cells

dn4<-susbet(dn4.int.reg, idents=c('DN4_1','DN4_2','DN4_3','DN4_4','DN4_5'))

DefaultAssay(dn4)<-'RNA'
Idents(dn4)<-dn4$geno

## Check Zfp335 expression
FindMarkers(dn4, ident.1 = 'Zfp335cKO', features = 'Zfp335', assay='RNA', logfc.threshold = 0)
#p_val  avg_log2FC pct.1 pct.2    p_val_adj
#Zfp335 1.017059e-34 -0.05482104  0.08 0.177 1.630651e-30

##Zfp335 expression not reliable for identifying true Zfp335 mutant cells
## Determine Zfp335 TF activity to identify mutant cells lacking Zfp335 using list of Zfp335 target genes down-regulated in Zfp335-deficient DP thymocytes

zfpTargetsDown<-read.RDS('Zfp335.Targets.DP.RNAseq.down.genes.rds')
ZfpTargetsDown<-intersect(ZfpTargetsDown,dn4@assays$RNA@counts@Dimnames[[1]])

ZfpTFscore<-100*(Matrix::colSums(dn4@assays$RNA@counts[ZfpTargetsDown,])/Matrix::colSums(dn4@assays$RNA@counts))

dn4<-AddMetaData(dn4,ZfpTFscore,'ZfpTargetsDown')
VlnPlot(dn4,'ZfpTargetsDown', group.by = 'geno', pt.size = 0, cols = c('grey50','red'))

### Determine cut off values for Zfp335 activity in mutant cells
library(cutoff)

meta<-dn4@meta.data
zmeta<- meta[meta$geno == 'Zfp335cKO',]

### zmeta$ZfpTargetsDown is gene signature score for Zfp335 target genes down-regulated in mutant DP thymocytes

length(zmeta$ZfpTargetsDown)

hist(zmeta$ZFPtargetsDown, 100, F, xlab = 'Zfp335 targets down-regulated in mutant DP thymocytes', ylab = 'density',main = NULL, col='grey')
lines(density(zmeta$ZFPtargetsDown), lwd=1.5, col='blue')

(td_out<-em(zmeta$ZFPtargetsDown,'normal','normal'))

#Finite mixture model fitting to dataset "zmeta$ZFPtargetsDown":

#  probability to belong to distribution 1 = 0.4524258 
#distribution 1: normal 
#location (mean) = 0.0893388
#scale (sd) = 0.0402361
#distribution 2: normal 
#location (mean) = 0.2588598
#scale (sd) = 0.042929
#deviance of the fitted model: -5158.4885047 

confint(td_out,nb=100,level = 0.95)

#         Estimate      2.5 %     97.5 %
#  mu1    0.08933881 0.08704260 0.09165061
# sigma1 0.04023607 0.03862047 0.04193928
# mu2    0.25885983 0.25663488 0.26107526
# sigma2 0.04292901 0.04137650 0.04456587
# lambda 0.45222814 0.43426205 0.47019423

hist(zmeta$ZFPtargetsDown, 100,F,xlab = 'Zfp335 targets down-regulated in mutant DP thymocytes', ylab = 'density',main = NULL, col='grey', ylim = c(0,10))
lines(td_out,lwd=1.5, col='red')

(td_cutoff <- cutoff(td_out))

#   Estimate     2.5 %    97.5 % 
#  0.2000175 0.1990926 0.2009425 

polygon(c(td_cutoff[-1],rev(td_cutoff[-1])),c(0,0,10,10), col=rgb(0,0,1,.2),border=NA)
abline(v=td_cutoff[-1], lty=2, col='blue')
abline(v=td_cutoff[1],col='blue')

### Use a cutoff value of 0.2000175 to ID cells with low Zfp335 TF activity

dn4[['genoTD']]<-paste(dn4$geno,ifelse(dn4$ZfpTargetsDown <= 0.2000175,'lo','hi'), sep = "_")

VlnPlot(dn4,'ZfpTargetsDown', group.by = 'genoTD')
Idents(dn4)<-dn4$genoTD

mut.markers<-FindMarkers(dn4,ident.1 = c('Zfp335cKO_lo','Zfp335cKO_hi'),ident.2 = c('WT_lo','WT_hi'), logfc.threshold = 0.05, assay = 'RNA')
mut.lo.markers<-FindMarkers(dn4,ident.1 = 'Zfp335cKO_lo',ident.2 = c('WT_lo','WT_hi'), logfc.threshold = 0.05, assay = 'RNA')
mut.hi.markers<-FindMarkers(dn4,ident.1 = 'Zfp335cKO_hi',ident.2 = c('WT_lo','WT_hi'), logfc.threshold = 0.05, assay = 'RNA')
mut.hi.vs.lo.markers<-FindMarkers(dn4,ident.1 = 'Zfp335cKO_lo', ident.2 = 'Zfp335cKO_hi', assay = 'RNA', logfc.threshold = 0.1)

EnhancedVolcano(mut.lo.markers,lab = rownames(mut.lo.markers), x='avg_log2FC', y='p_val_adj', pCutoff = 1e-50,FCcutoff = 0.15, col = c('grey50','grey50','grey50','red'), xlim = c(-0.75,0.75))
EnhancedVolcano(mut.hi.markers,lab = rownames(mut.hi.markers), x='avg_log2FC', y='p_val_adj', pCutoff = 1e-50,FCcutoff = 0.15, col = c('grey50','grey50','grey50','red'), xlim=c(-0.7,0.7))
EnhancedVolcano(mut.hi.vs.lo.markers,lab = rownames(mut.hi.vs.lo.markers), x='avg_log2FC', y='p_val_adj', pCutoff = 1e-50,FCcutoff = 0.15, col = c('grey50','grey50','grey50','red'), xlim=c(-0.7,0.7))


mut.lo.filtered<-filter(mut.lo.markers, abs(avg_log2FC) >= 0.15 & p_val_adj <= 1e-50)
mut.hi.filtered<-filter(mut.hi.markers, abs(avg_log2FC) >= 0.15 & p_val_adj <= 1e-50)
mut.hi.vs.lo.markers.filtered<-filter(mut.hi.vs.lo.markers, abs(avg_log2FC) >= 0.15 & p_val_adj <= 1e-50)

### Subset data to remove non-mutant Zfp335cKO cells

truMutDN4<- subset(dn4, idents=c('WT_hi','WT_lo','Zfp335cKO_lo'))
truMutDN4<- ScaleData(truMutDN4, assay='RNA', model.use = 'negbinom')
truMutDN4<-FindVariableFeatures(truMutDN4, selection.method = 'vst',nfeatures = 5000, assay = 'RNA')
length(truMutDN4@assays$RNA@var.features)

### Impute dropouts using ALRA
truMutDN4<- RunALRA(truMutDN4, assay='RNA', features=truMutDN4@assays$RNA@counts@Dimnames[[1]])

### Recluster Cells
truMutDN4<-ScaleData(truMutDN4)
truMutDN4<-RunPCA(truMutDN4,features=truMutDN4@assays$RNA@var.features, nfeatures.print=10)


dn4<-RunPCA(dn4,features=VariableFeatures(dn4),nfeatures.print = 10)
truMutDN4<-JackStraw(truMutDN4,num.replicate = 100,dims = 40)
truMutDN4<-ScoreJackStraw(truMutDN4,dims=1:40)
JackStrawPlot(truMutDN4, dims=1:40)

truMutDN4<-FindNeighbors(truMutDN4, dims=1:40)
truMutDN4<-FindClusters(truMutDN4, resolution = 0.5)
truMutDN4<-RunUMAP(truMutDN4, dims=1:40)
DimPlot(truMutDN4, group.by = 'geno',pt.size = 1.5, cols = c('grey50','red'))
DimPlot(truMutDN4,split.by = 'geno', pt.size = 1, label = T, label.size = 5)
DimPlot(truMutDN4,group.by = 'Phase', pt.size = 1)

table(truMutDN4$seurat_clusters, truMutDN4$geno)

#    WT Zfp335cKO
#0  876       419
#1 1149       103
#2  824        37
#3  672        52
#4    1       722
#5  150        65
#6   29        77
#7   67        34

table(truMutDN4$Phase)
#G1  G2M    S
#0    0 1015  280
#1   11  289  952
#2    2    2  857
#3    1  147  576
#4    2   86  635
#5   16   53  146
#6    6   47   53
#7   84    1   16

table(truMutDN4$geno)

#   WT Zfp335cKO 
#   3768      1509

table(truMutDN4$seurat_clusters,truMutDN4$geno)

#WT Zfp335cKO
#0  876       419
#1 1149       103
#2  824        37
#3  672        52
#4    1       722
#5  150        65
#6   29        77
#7   67        34

## Run differential expression on each cluster
dn4.markers<-FindAllMarkers(truMutDN4,assay = 'RNA')

dn4.markers %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_logFC) -> top25

# Export list of top 25 marker genes from each cluster for GO analysis using DAVID

write.csv(top25,'TrueMutant.DN4.top25.markers.by.cluster.csv')

# DAVID analysis showed nearly all enriched pathways among the top cluster-defining genes were realted to cell cycle

## To get better understanding of Zfp335 regulated genes in dataset compare all WT to all Zfp335cKO cells
Idents(truMutDN4)<-truMutDN4$geno
LoMutMarkers<-FindMarkers(truMutDN4,ident.1 = 'Zfp335cKO', ident.2 = 'WT', logfc.threshold = 0.05, min.pct = 0.1, assay = 'RNA')

## Identify differentially expressed Zfp335 target genes in Zfp335cKO_lo cells
## List of Zfp335 target genes determined from previous Zfp335 thymus ChIP-seq data (GSE)

thy<-c("1110038B12Rik",	"1700034F02Rik",	"1700063D05Rik",	"1700071M16Rik",	"2410002F23Rik",	"4932411N23Rik",	"A530054K11Rik",	"AA987161",	"Actr10",	"Aimp1",
       "Alg2",	"Alkbh5",	"Anapc1",	"Anapc7",	"Ankle2",	"Ap3m1",	"Arhgap11a",	"Atp5l",	"Bag6",	"BC018507",
       "Bckdha",	"Birc2",	"Blcap",	"Bud31",	"C130046K22Rik",	"C920009B18Rik",	"Cactin",	"Caprin1",	"Ccdc43",	"Ccdc58",
       "Cct4",	"Cdc73",	"Cep76",	"Chd8",	"Cnpy2",	"Cox7a2l",	"Cs",	"D11Wsu47e",	"Dand5",	"Dctn1",
       "Ddx31",	"Ddx39b",	"Dgkz",	"Dnal4",	"Dph6",	"Dync1li2",	"Efnb1",	"Eftud2",	"Eif2s3x",	"Eif4enif1",
       "Eif4h",	"Fam222b",	"Fars2",	"Farsb",	"Fosb",	"Gm21119",	"Gm21276",	"Gorab",	"H2-M10.1",	"Hadha",
       "Hilpda",	"Hspa9",	"Ipo13",	"Iqcg",	"Irgq",	"LOC100503496",	"Lrr1",	"Mapkapk5",	"Mir592",	"Mkks",
       "Mrpl11",	"Mrpl38",	"Mrpl9",	"Mrps5",	"Mrps7",	"Mrto4",	"Mtss1",	"Nat10",	"Ncln",	"Ndc1",
       "Ndufa10",	"Ndufa4",	"Nek10",	"Nhlrc3",	"Nme6",	"Nsun2",	"Nudcd3",	"Nup205",	"Nup85",	"Nup98",
       "Nxph1",	"Olfr161",	"Pcbd2",	"Pcbp1",	"Pdap1",	"Pes1",	"Pi4kb",	"Pin4",	"Pno1",	"Polr2e",
       "Polr3f",	"Psmd1",	"Psmd12",	"Psmg2",	"Ptges2",	"Rabggtb",	"Rabl2",	"Rad17",	"Rbbp5",	"Rbks",
       "Rexo1",	"Rpl21",	"Rpl35a",	"Rpl7",	"Rpn1",	"Rps10",	"Rps12",	"Rps23",	"Rps24",	"Rps27a",
       "Sart3",	"Sec11c",	"Selrc1",	"Sep15",	"Sf3b3",	"Skp1a",	"Slc35b1",	"Smc2",	"Spcs1",	"Stam2",
       "Stamos",	"Suclg1",	"Suds3",	"Supv3l1",	"Taf9",	"Tm9sf1",	"Tmed1",	"Tmppe",	"Tomm20",	"Trap1",
       "Tsn",	"Ttc8",	"Tubd1",	"Txnrd2",	"Ube2i",	"Ufm1",	"Upf2",	"Urod",	"Usp1",	"Ust",
       "Vars",	"Vdac2",	"Vti1b",	"Wdr47",	"Wdr92",	"Xpnpep3",	"Zfp335",	"Zfp383",	"Zfp541",	"Zfp738",
       "Zranb2",	"Hadhb")
thy<-intersect(thy, truMutDN4@assays$RNA@counts@Dimnames[[1]])

LoMutZfpTargets<-na.omit(LoMutMarkers[thy,])
EnhancedVolcano(LoMutZfpTargets,lab=rownames(LoMutZfpTargets),'avg_log2FC','p_val_adj', FCcutoff = 0.1, xlim = c(-1,0.25), col=c('grey50','grey50','grey50','red'))
SigLoMutZfpTargets<-filter(LoMutZfpTargets, abs(avg_log2FC) >=0.15 & p_val_adj < 1e-50)

## Identify Zfp335 targets involved with regulating cell death
### Downloaded list of all mouse genes with experimental evidence of negative regulation in apoptotic processes from: https://www.ebi.ac.uk/QuickGO/term/GO:0043066
### Selected Mus musculus genes and filtered by annotations wit evidence by manual assertion
### Filtered output in excel to remove duplicates and create gene list


NegApop<-readRDS("NegativeRegulatorsOfApoptosis.Mouse.RDS")
NegApop<-intersect(NegApop,truMutDN4@assays$RNA@counts@Dimnames[[1]])
ZfpNegApop<-intersect(NegApop,thy)

ZfpNegApop.Markers<-FindMarkers(truMutDN4,ident.1 = 'Zfp335cKO',ident.2 = 'WT', assay = 'RNA',features = ZfpNegApop,logfc.threshold = 0, min.pct = 0.01)
EnhancedVolcano(ZfpNegApop.Markers,lab = rownames(ZfpNegApop.Markers),x='avg_log2FC',y='p_val_adj',FCcutoff = 0.15, pCutoff = 1e-50,xlim = c(-0.8,0.1),drawConnectors = T, col=c('grey50','grey50','grey50','red'), selectLab = ZfpNegApop)
ZfpNegApop.Markers.ALRA<-FindMarkers(truMutDN4,ident.1 = 'Zfp335cKO',ident.2 = 'WT', assay = 'alra',features = ZfpNegApop,logfc.threshold = 0, min.pct = 0.01)
ZfpNegApop.Markers.ALRA[['dPct']]<-ZfpNegApop.Markers.ALRA$pct.1 - ZfpNegApop.Markers.ALRA$pct.2


#### Start unbiased pathway analysis using SingleCellSignatureExplorer to compute cell-by-cell geneset scores using the reactome pathway database
library(nichenetr)
library(tidyverse)
counts<-t(as.matrix(truMutDN4@assays$RNA@data))
colnames(counts)<- counts %>% colnames() %>% convert_mouse_to_human_symbols()
counts<- counts %>% .[!is.na(rownames(counts)), !is.na(colnames(counts))]
head(colnames(counts))
write.table(counts,'tab.tsv',sep = "\t", quote = F)

## Downloaded the reactome database genesets from MSigDb and performed geneset scoring using SingleCellSignatureScorer
## Read in geneset scores matrix output by SingleCellSignatureScorer

scores<-read.csv('C2_CP_REACTOME_tab.tsv',sep="\t", row.names = 1)

## Add scores to seurat object as an assay

truMutDN4[['reactome']]<-CreateAssayObject(data=t(scores))

## Run wilcoxon tests to compare Mut to WT genesets

reactome.markers<-FindMarkers(truMutDN4, ident.1 = 'Zfp335cKO', ident.2 = 'WT', assay='reactome', logfc.threshold = 0)

## Make viridis color scale for split FeaturePlots
vir<-c("#440154FF","#481567FF","#482677FF","#453781FF","#404788FF","#39568CFF","#33638DFF","#2D708EFF","#238A8DFF","#1F968BFF","#20A387FF","#29AF7FFF","#3CBB75FF","#55C667FF","#73D055FF","#95D840FF","#B8DE29FF","#DCE319FF","#FDE725FF")

FeaturePlot(truMutDN4,'REACTOME-IRF3-MEDIATED-ACTIVATION-OF-TYPE-1-IFN', pt.size = 1.5, split.by = 'geno', max.cutoff = 'q95', cols = vir)

### Make violin plots for STING-activated pro-apoptotic Bcl2 genes, recolored violins in inkscape

VlnPlot(truMutDN4,c('Bbc3','Pmaip1','Bcl2l11','Bax'), pt.size = 0, stack = T)+NoLegend()

##### Additional analyses for revised manuscript

## Make new violin plot for Fig 4E 

VlnPlot(dn4.int.reg, c('Nr4a1','Cd69','Pdcd1','Egr1','Cd2','Itm2a','Cd24a'), pt.size=0, stack=T, flip=T, group.by='cellType', fill.by='ident')+NoLegend()

#### Make violin plots of IFN-I signaling and mTOR singaling gene signatures for 'true' mutant DN4 data

VlnPlot(truMutDN4, c('REACTOME-INTERFERON-ALPHA-BETA-SIGNALING','REACTOME-MTOR-SIGNALLING'), pt.size=0, group.by='seurat_clusters', split.by='geno', split.plot=T, cols=c('grey50','red'))

#### Compare signature scores between genotypes within each cluster

truMutDN4[['genoCT']]<-paste(truMutDN4$geno, truMutDN4$seurat_clusters, sep='_')
Idents(truMutDN4)<-truMutDN4$genoCT

sigs<-c('REACTOME-INTERFERON-ALPHA-BETA-SIGNALING','REACTOME-MTOR-SIGNALLING')

FindMarkers(truMutDN4, ident.1='Zfp335cKO_0', ident.2='WT_0', features=sigs, assay='reactome, logfc.threshold=0)
#                                                p_val  avg_log2FC pct.1 pct.2  p_val_adj
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 1.313689e-05 0.018347257     1     1 0.01998121
#REACTOME-MTOR-SIGNALLING                 4.158618e-02 0.007083467     1     1 1.00000000

FindMarkers(truMutDN4, ident.1='Zfp335cKO_1', ident.2='WT_1', features=sigs, assay='reactome, logfc.threshold=0)
#                                               p_val    avg_log2FC pct.1 pct.2 p_val_adj
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 0.000982289  0.0213752721     1     1         1
#REACTOME-MTOR-SIGNALLING                 0.999092070 -0.0002743283     1     1         1

FindMarkers(truMutDN4, ident.1='Zfp335cKO_2', ident.2='WT_2', features=sigs, assay='reactome, logfc.threshold=0)
#                                              p_val avg_log2FC pct.1 pct.2 p_val_adj
#REACTOME-MTOR-SIGNALLING                 0.07414677 0.01237910     1     1         1
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 0.10644527 0.01666793     1     1         1

FindMarkers(truMutDN4, ident.1='Zfp335cKO_3', ident.2='WT_3', features=sigs, assay='reactome, logfc.threshold=0)
#                                             p_val  avg_log2FC pct.1 pct.2 p_val_adj
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 0.2200214  0.01437493     1     1         1
#REACTOME-MTOR-SIGNALLING                 0.5617761 -0.00906888     1     1         1

FindMarkers(truMutDN4, ident.1='Zfp335cKO_4', ident.2='WT_4', features=sigs, assay='reactome, logfc.threshold=0)
#Error in ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2,  : 
#  Cell group 2 has fewer than 3 cells

FindMarkers(truMutDN4, ident.1='Zfp335cKO_5', ident.2='WT_5', features=sigs, assay='reactome, logfc.threshold=0)
#                                              p_val   avg_log2FC pct.1 pct.2 p_val_adj
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 0.01422043 0.0334002745     1     1         1
#REACTOME-MTOR-SIGNALLING                 0.40207696 0.0007438672     1     1         1

FindMarkers(truMutDN4, ident.1='Zfp335cKO_6', ident.2='WT_6', features=sigs, assay='reactome, logfc.threshold=0)
#                                             p_val avg_log2FC pct.1 pct.2 p_val_adj
#REACTOME-MTOR-SIGNALLING                 0.1815682 0.04857284     1     1         1
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 0.3477086 0.03990863     1     1         1

FindMarkers(truMutDN4, ident.1='Zfp335cKO_7', ident.2='WT_7', features=sigs, assay='reactome, logfc.threshold=0)
#                                              p_val  avg_log2FC pct.1 pct.2 p_val_adj
#REACTOME-INTERFERON-ALPHA-BETA-SIGNALING 0.03617595  0.04672645     1     1         1
#REACTOME-MTOR-SIGNALLING                 0.81810777 -0.00805391     1     1         1


