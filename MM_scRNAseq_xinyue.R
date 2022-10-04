fun_packages <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      #Si le package ne peut pas être chargé, on l'installe
      install.packages( i , dependencies = TRUE )
      #Charger le package après l'installation
      require( i , character.only = TRUE )
    }
    if (!requireNamespace("BiocManager", quietly = TRUE))
    {install.packages("BiocManager")}
  }
}

fun_packages(c("Seurat","stringr","ggplot2","SingleR","celldex"))

#BiocManager::install("Seurat")
#BiocManager::install("celldex")
#BiocManager::install("SingleR")

#load data of scRNAseq
fs=list.files('./GSM_file/')

GSM_files_list=lapply(fs, function(x){
  a=read.table(file.path('./GSM_file/',x))
})

GSM_files_count=do.call(cbind,GSM_files_list)
GSM_files_count=as.data.frame(GSM_files_count)
#load the metadata which contains the infos of scRNAseq
metadata=read.table(gzfile("GSE161195_metadata_s.txt.gz"))

#effet batch
#list_test=list()
#test=function(GSM_files_list){
#for (i in 1:length(GSM_files_list)){
#colnames(GSM_files_list[[i]])=gsub("_", "-", colnames(GSM_files_list[[i]]))
#rownames(GSM_files_list[[i]])=gsub("_", "-", rownames(GSM_files_list[[i]]))
#seuratObj_test=CreateSeuratObject(counts = GSM_files_list[[i]],meta.data=metadata)
#seuratObj_test[["percent.mt"]] <- PercentageFeatureSet(seuratObj_test, pattern = "^MT-")
#seuratObj_subset <- subset(seuratObj_test, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 50)
#seuratObj_subset <- NormalizeData(seuratObj_subset, normalization.method = "LogNormalize", scale.factor = 10000)
#seuratObj_subset_variable <- FindVariableFeatures(seuratObj_subset, selection.method = "vst", nfeatures = 2000)
#list_test[[i]]=seuratObj_subset_variable
#}
#return (list_test)
#}
#list_yyy=test(GSM_files_list)

#mmrf <- FindIntegrationAnchors(object.list = list_yyy, dims = 1:10)
#mmrf.integrated <- IntegrateData(anchorset = mmrf, dims = 1:10)

#preprocessing of data
colnames(GSM_files_count)=gsub("_", "-", colnames(GSM_files_count))
rownames(GSM_files_count)=gsub("_", "-", rownames(GSM_files_count))
#construction of data object
seuratObj_test=CreateSeuratObject(counts = GSM_files_count, project = "MM",meta.data=metadata)

#observation of pourcentage of mitochondria and count of RNA so we can filter with sandard
seuratObj_test[["percent.mt"]] <- PercentageFeatureSet(seuratObj_test, pattern = "^MT-")
head(seuratObj_test@meta.data, 5)
VlnPlot(seuratObj_test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seuratObj_test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seuratObj_test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filter of data and observation of distrubution of data
seuratObj_subset <- subset(seuratObj_test, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 50)
hist(colSums(seuratObj_subset$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

#normalization of data
seuratObj_subset <- NormalizeData(seuratObj_subset, normalization.method = "LogNormalize", scale.factor = 10000)
hist(colSums(seuratObj_subset$RNA@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")  

seuratObj_subset_variable <- FindVariableFeatures(seuratObj_subset, selection.method = "vst", nfeatures = 2000)


# mark top 10 expressed genes (just important for this data)
top10 <- head(VariableFeatures(seuratObj_subset_variable), 10)
#[1] "HBB"     "HBA2"    "IGHM"    "HBA1"    "CXCL12"  "IGHD"   
#[7] "MIR4539" "IGLC7"   "IGHG4"   "IGLC6" 
all.genes <- rownames(seuratObj_subset_variable)
#scale the data
seuratObj_subset_scale <- ScaleData(seuratObj_subset_variable, features = all.genes)
#pca
seuratObj_PCA <- RunPCA(seuratObj_subset_scale, features = VariableFeatures(object = seuratObj_subset_scale))
#the calculate result of pca
print(seuratObj_PCA[["pca"]], dims = 1:5, nfeatures = 5)

#observer the best dimension of pca
ElbowPlot(object = seuratObj_PCA,
          ndims = 40)

#calculate for different clusters
#method of Umap
seuratObj_cluster <- seuratObj_PCA %>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunUMAP(dims = 1:15) %>% 
  identity()

#method of tsne
seuratObj_cluster <- seuratObj_PCA %>% 
  RunTSNE(reduction = "pca", dims = 1:15)%>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = 0.5)

#visualization of data cluster
DimPlot(seuratObj_cluster, reduction = "umap",label=T)
DimPlot(seuratObj_cluster, reduction = "tsne",label=T)

############# Annotation cluster with SingleR ###########
DefaultAssay(seuratObj_cluster) <- "RNA"
# for B cells :  cluster, 1,21

sce_for_SingleR <- GetAssayData(seuratObj_cluster, slot="data")
clusters=seuratObj_cluster@meta.data$seurat_clusters
Immune.se=DatabaseImmuneCellExpressionData()

pred.mmrf_immuno <- SingleR(test = sce_for_SingleR, ref = Immune.se, 
                            labels = Immune.se$label.main)

table(pred.mmrf_immuno$labels,seuratObj_cluster@meta.data$seurat_clusters)
#             
#                0     1     2     3     4     5     6     7     8     9    10
#B cells       12881  8083  3897  3064  2654  2000  2258    77  1552  1462   220
#Monocytes       428    13   350   724    33   521     4  2019     0     2   180
#NK cells        109     1   134   241     3    44     4     4     1     0   278
#T cells, CD4+    55     2    76   286     2    35     0     2     0     0   641
#T cells, CD8+    34     0    58   121     0    17     1     2     0     0    73

#                11    12    13    14    15    16
#B cells        1139   593   963   905   953    10
#Monocytes         0   366    22    62     8   436
#NK cells          0    15     8     6     3     6
#T cells, CD4+     0    46     3     3     9     0
#T cells, CD8+     0    17     6     9     1     0

#by seeing this table, we determmine the nature of data

cluster.ids <- c("B cells", "B cells", "B cells", "B cells", "B cells", "B cells","B cells",
                 "Monocytes", "B cells", "B cells","T cells, CD4+","B cells",
                 "B cells","B cells","B cells","B cells","Monocytes")

#annotate the cluster according to this result
names(cluster.ids) <- levels(seuratObj_cluster)
seuratObj_cluster <- RenameIdents(seuratObj_cluster, cluster.ids)
DimPlot(seuratObj_cluster, reduction = "umap", label = TRUE, pt.size = 0.5)

#the automate annotation without the distinction
#seuratObj_cluster@meta.data$labels=pred.mmrf_immuno$labels
#DimPlot(seuratObj_cluster, group.by ="labels",reduction = "umap",label=T)

gene_list=c("SNORD89",  "AURKAIP1", "RNVU1-18","ELANE", "LINC02108","ANKRD11",
            "NCOA3","CYLD","RNY1","RFESD","RNU1-1","RNU1-3")

gene_biblio=c ("KRAS","NRAS","FAM46C","DIS3","TP53", "BRAF", "TRAF3", "PRDM1", "CYLD", 
               "RB1", "IRF4", "EGR1","HIST1H1E","ACTG1")

#VlnPlot(object = seuratObj_cluster, features =gene_biblio,log =T ) 
FeaturePlot(seuratObj_cluster, features = gene_biblio)
FeaturePlot(seuratObj_cluster, features = gene_list)

#regarding the plot generated, we select finally
#AURKAIP1,ANKRD11, NCOA3, CYLD, FAM46C,PRDM1,IRF4,ACTG1,
#as our potentiel biomarkers for survival analysis
