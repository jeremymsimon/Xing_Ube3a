library(reticulate)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(sctransform)
library(fishpond)
library(tximport)
library(devtools)

# path to the Alevin output files
allfiles = c("e14_T485A1_1_SM0104_N714_alevin_quants_mat.gz","e14_T485A1_SM0104_N705_alevin_quants_mat.gz","e14_T485A2_SM0104_N706_alevin_quants_mat.gz","e14_T485A3_SM0104_N707_alevin_quants_mat.gz","e14_T485A4_SM0104_N710_alevin_quants_mat.gz","e14_WT10_SM0104_N703_alevin_quants_mat.gz","e14_WT11_SM0104_N704_alevin_quants_mat.gz","e14_WT8_1_SM0104_N711_alevin_quants_mat.gz","e14_WT8_2_SM0104_N701_alevin_quants_mat.gz","e14_WT9_1_SM0104_N712_alevin_quants_mat.gz","e14_WT9_2_SM0104_N702_alevin_quants_mat.gz","p0_T485A1_SM0104_N719_alevin_quants_mat.gz","p0_T485A2_SM0104_N720_alevin_quants_mat.gz","p0_T485A3_SM0104_N721_alevin_quants_mat.gz","p0_T485A4_SM0104_N722_alevin_quants_mat.gz","p0_T485A5_SM0104_N723_alevin_quants_mat.gz","p0_WT2_SM0104_N715_alevin_quants_mat.gz","p0_WT3_SM0104_N716_alevin_quants_mat.gz","p0_WT4_SM0104_N718_alevin_quants_mat.gz")
allsamples = c("e14_T485A1_1","e14_T485A1","e14_T485A2","e14_T485A3","e14_T485A4","e14_WT10","e14_WT11","e14_WT8_1","e14_WT8_2","e14_WT9_1","e14_WT9_2","p0_T485A1","p0_T485A2","p0_T485A3","p0_T485A4","p0_T485A5","p0_WT2","p0_WT3","p0_WT4")
allsamples=as.list(allsamples)

#Create list of samples
mylist = allsamples
names(mylist) = unlist(mylist)

#Create objects for each sample
for(i in 1:length(allfiles)){
	file <- file.path(allfiles[i])
	txi <- tximport(file, type="alevin")
	txi.annot = txi$counts
	name = gsub("_","-",allsamples[[i]])
	colnames(txi.annot) = paste0(name,"_",colnames(txi.annot))
	assign(allsamples[[i]], CreateSeuratObject(counts = txi.annot , min.cells = 3, min.features = 300, project = allsamples[[i]]))
}

#Create list of sample objects
ctx.all <- merge(e14_T485A1_1, y = c(e14_T485A1,e14_T485A2,e14_T485A3,e14_T485A4,e14_WT10,e14_WT11,e14_WT8_1,e14_WT8_2,e14_WT9_1,e14_WT9_2,p0_T485A1,p0_T485A2,p0_T485A3,p0_T485A4,p0_T485A5,p0_WT2,p0_WT3,p0_WT4), project = "e14.p0.wt.t485a.ctx")
ctx.list <- SplitObject(object = ctx.all, split.by = "ident")


#Add Replicate Metadata
p0_WT2@meta.data$replicate <- "1"
p0_WT3@meta.data$replicate <- "2"
p0_WT4@meta.data$replicate <- "3"
p0_T485A1@meta.data$replicate <- "1"
p0_T485A2@meta.data$replicate <- "2"
p0_T485A3@meta.data$replicate <- "3"
p0_T485A4@meta.data$replicate <- "4"
p0_T485A5@meta.data$replicate <- "5"
e14_WT10@meta.data$replicate <- "1"
e14_WT11@meta.data$replicate <- "2"
e14_WT8_1@meta.data$replicate <- "3"
e14_WT8_2@meta.data$replicate <- "4"
e14_WT9_1@meta.data$replicate <- "5"
e14_WT9_2@meta.data$replicate <- "6"
e14_T485A1_1@meta.data$replicate <- "1"
e14_T485A1@meta.data$replicate <- "2"
e14_T485A2@meta.data$replicate <- "3"
e14_T485A3@meta.data$replicate <- "4"
e14_T485A4@meta.data$replicate <- "5"

#Add genotype metadata
p0_WT2@meta.data$gt <- "WT"
p0_WT3@meta.data$gt <- "WT"
p0_WT4@meta.data$gt <- "WT"
p0_T485A1@meta.data$gt <- "T485A"
p0_T485A2@meta.data$gt <- "T485A"
p0_T485A3@meta.data$gt <- "T485A"
p0_T485A4@meta.data$gt <- "T485A"
p0_T485A5@meta.data$gt <- "T485A"
e14_WT10@meta.data$gt <- "WT"
e14_WT11@meta.data$gt <- "WT"
e14_WT8_1@meta.data$gt <- "WT"
e14_WT8_2@meta.data$gt <- "WT"
e14_WT9_1@meta.data$gt <- "WT"
e14_WT9_2@meta.data$gt <- "WT"
e14_T485A1_1@meta.data$gt <- "T485A"
e14_T485A1@meta.data$gt <- "T485A"
e14_T485A2@meta.data$gt <- "T485A"
e14_T485A3@meta.data$gt <- "T485A"
e14_T485A4@meta.data$gt <- "T485A"

#Add timepoint metadata
p0_WT2@meta.data$time <- "P0"
p0_WT3@meta.data$time <- "P0"
p0_WT4@meta.data$time <- "P0"
p0_T485A1@meta.data$time <- "P0"
p0_T485A2@meta.data$time <- "P0"
p0_T485A3@meta.data$time <- "P0"
p0_T485A4@meta.data$time <- "P0"
p0_T485A5@meta.data$time <- "P0"
e14_WT10@meta.data$time <- "E14"
e14_WT11@meta.data$time <- "E14"
e14_WT8_1@meta.data$time <- "E14"
e14_WT8_2@meta.data$time <- "E14"
e14_WT9_1@meta.data$time <- "E14"
e14_WT9_2@meta.data$time <- "E14"
e14_T485A1_1@meta.data$time <- "E14"
e14_T485A1@meta.data$time <- "E14"
e14_T485A2@meta.data$time <- "E14"
e14_T485A3@meta.data$time <- "E14"
e14_T485A4@meta.data$time <- "E14"



#Flter and normalize all samples
for (i in 1:length(x = ctx.list)) {
	ctx.list[[i]] <- PercentageFeatureSet(object = ctx.list[[i]], pattern = "^HUMAN-MT-", col.name = "percent.hum.mt")
	ctx.list[[i]] <- PercentageFeatureSet(object = ctx.list[[i]], pattern = "^mt-", col.name = "percent.mt")
    ctx.list[[i]] <- PercentageFeatureSet(object = ctx.list[[i]], pattern = "^Hbb-", col.name = "percent.Hbb")
	ctx.list[[i]] <- PercentageFeatureSet(object = ctx.list[[i]], pattern = "^HUMAN-", col.name = "percent.human")
    ctx.list[[i]] <- subset(ctx.list[[i]], subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5 & percent.Hbb < 5 & percent.hum.mt < 5 & percent.human < 5)
	ctx.list[[i]] <- SCTransform(ctx.list[[i]], vars.to.regress = c("percent.mt","percent.hum.mt","percent.Hbb"), verbose = FALSE)
}



# Set up integration, then integrate data
### NOTE: I turned verbose to TRUE so I can monitor the process
options(future.globals.maxSize=5242880000)
ctx.features <- SelectIntegrationFeatures(object.list = ctx.list, nfeatures = 5000)
ctx.list <- PrepSCTIntegration(object.list = ctx.list, anchor.features = ctx.features, verbose = TRUE)
ctx.anchors <- FindIntegrationAnchors(object.list = ctx.list, normalization.method = "SCT", anchor.features = ctx.features, verbose = TRUE)
ctx.integrated <- IntegrateData(anchorset = ctx.anchors, normalization.method = "SCT", verbose = TRUE)


# Reset metadata
ctx.integrated@meta.data$time = toupper(unlist(lapply(strsplit(ctx.integrated@meta.data$orig.ident, "-"), '[[', 1)))
ctx.integrated@meta.data$gt = gsub("\\d+$","",toupper(unlist(lapply(strsplit(ctx.integrated@meta.data$orig.ident, "-"), '[[', 2))),perl=T)
repInfo=c(1,2,3,1,2,3,4,5,1,2,3,4,5,6,1,2,3,4,5)
names(repInfo) = c("p0_WT2","p0_WT3","p0_WT4","p0_T485A1","p0_T485A2","p0_T485A3","p0_T485A4","p0_T485A5","e14_WT10","e14_WT11","e14_WT8_1","e14_WT8_2","e14_WT9_1","e14_WT9_2","e14_T485A1_1","e14_T485A1","e14_T485A2","e14_T485A3","e14_T485A4")

for(i in 1:length(ctx.integrated@meta.data$orig.ident)) {
	name = gsub("-","_",ctx.integrated@meta.data$orig.ident[i])
	rep = as.numeric(repInfo[name])
	ctx.integrated@meta.data$replicate[i] = rep
}


#Run PCA
ctx.integrated <- RunPCA(ctx.integrated, verbose = FALSE, npcs = 100)
ctx.integrated <- RunUMAP(ctx.integrated, dims = 1:100,umap.method="umap-learn",metric="correlation")
ctx.integrated <- FindNeighbors(ctx.integrated, dims = 1:100, verbose = FALSE)
ctx.integrated <- FindClusters(ctx.integrated, verbose = FALSE, resolution = 2, algorithm=2)


#Plot UMAPs
DimPlot(ctx.integrated,reduction = "umap", label = TRUE)
DimPlot(ctx.integrated,reduction = "umap", group.by = "gt")
DimPlot(ctx.integrated,reduction = "umap", group.by = "time")
DimPlot(ctx.integrated,reduction = "umap", group.by = "replicate")


DimPlot(ctx.integrated,reduction = "umap", label = TRUE, split.by = "gt")
DimPlot(ctx.integrated,reduction = "umap", label = TRUE, split.by = "time")
DimPlot(ctx.integrated,reduction = "umap", label = TRUE, split.by = "replicate")





# QC plots

features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.hum.mt","percent.human")
plotwidth = 50
plotheight = 7

for (i in features) {
    print(VlnPlot(ctx.integrated, features = i, split.by = "gt"))
}





#Detect marker genes of each cluster
for (i in levels(ctx.integrated)) {
    markers <- FindMarkers(ctx.integrated, assay="SCT", slot="scale.data", ident.1=i)
    write.table(markers,paste0("scRNA_T485A_Markers_SCT_",i,".txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
}





#Make dotplots
goi <- rev(unique(c("Neurod6","Vim","Reln","Satb2","Rprm","Cntn2","Bcl11b","Fezf2","Mc4r","Crym","Nln","Lmo1","9130024F11Rik","Tfap2d","Cbln2","Nr2f1","Nr2f2","Lmo3","Sstr2","Eomes","Hes5","Sfrp1","Mki67","Ccnb1","Mcm3","Ung","Hmmr","Htra1","Gad1","Gad2","Dlx1","Dlx2","Lhx6","Sst","Zcchc12","Zic1","Pnoc","Adarb2","Cdca7","Syt6","Isl1","Six3","Tcf7l2","Slc1a3","Aldh1l1","Pax6","Olig2","Otx2","Trem2","Igfbp7","Egfl7","Igf2")))
DotPlot(ctx.integrated,features=goi,assay="SCT",cols = c("blue","red")) + RotatedAxis()



