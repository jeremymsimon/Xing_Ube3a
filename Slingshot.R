# Psuedotiming analysis with slingshot on interneuronal lineage
# Requires update of Seurat object to new version
# Now in R v4.1

library(umap)
library(princurve)
library(slingshot)
library(tidyverse)


subset = subset(ctx.integrated,subset = seurat_clusters %in% c(5,7,8,13,20,24,25,29,33))
subset = UpdateSeuratObject(object = subset)

# Find variably expressed/marker genes just of these cells
markers = FindAllMarkers(subset,assay="RNA",only.pos=T,logfc.threshold=1)

all_markers = unique(markers$gene)

# Create filtered data object on just this subset of cells and variable genes
subset.filtered = GetAssayData(subset,assay="integrated")
subset.filtered = subset.filtered[rownames(subset.filtered) %in% all_markers,]

# Re-compute UMAP reduction just on these cells and genes
subset.filtered.umap = umap::umap(t(as.matrix(subset.filtered)))
subset.umap = subset.filtered.umap$layout
colnames(subset.umap) = c("Dim1","Dim2")

# Run slingshot
subset.sling = slingshot(subset.umap, subset$seurat_clusters,end.clus=c(5,7))

# Retrieve clusters in lingeage order with
slingLineages(subset.sling)

# $Lineage1
# [1] "24" "29" "8"  "33" "20" "13" "25" "5" 
# 
# $Lineage2
#[1] "24" "29" "8"  "33" "20" "13" "7" 

# Plot using ggplot
# Following vignette here: http://www.bioconductor.org/packages/release/bioc/vignettes/traviz/inst/doc/slingshot.html

df = as.data.frame(cbind(subset.umap,"cl" = as.numeric(as.character(subset$seurat_clusters))))
curves <- slingCurves(subset.sling, as.df = TRUE)

ggplot(df, aes(x = -Dim2, y = Dim1)) +
  	geom_point(aes(color = factor(cl))) + 
  	theme_classic() +
	geom_path(data = curves %>% arrange(Order), aes(group = Lineage), size = 1, arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches"))) +
	guides(color=guide_legend(title="Cluster")) +
	xlab("UMAP 2 (reversed)") +
	ylab("UMAP 1") +
	scale_color_manual(breaks=c(24,29,8,33,20,13,7,25,5),values=c("#00B9E3", "#DB72FB", "#93AA00", "#FF61C3", "#00C19F", "#00BA38", "#D39200", "#619CFF", "#F8766D"),labels=c("24_GE_M", "29_GE_M", "08_GE_G1S", "33_GE_G1S", "20_Int4", "13_Int3", "07_Int1", "25_Int2_Zic+", "05_Int2"))




# Now prep for heatmap figure

library(tradeSeq)
library(RColorBrewer)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(presto)
library(viridis)

set.seed(123)

subset.gam = subset(ctx.integrated,subset = seurat_clusters %in% c(5,7,8,13,20,24,25,29,33), features=all_markers)
subset.gam = UpdateSeuratObject(object = subset.gam)

sce = as.SingleCellExperiment(subset.gam,assay = "RNA")

sce <- fitGAM(counts = counts(sce), sds = subset.sling)

ATres <- associationTest(sce, lineages = TRUE)

# Focus on Lineage 1, as this one includes Int2
topgenes <- rownames(ATres[order(ATres$pvalue_1), ])[1:250]


# Explore top marker genes  significantly temporally expressed in either lineage
ts.topgenes <- rownames_to_column(ATres,var="Gene") %>% 
	as_tibble() %>%
	filter(p.adjust(pvalue_1)<0.05 | p.adjust(pvalue_2)<0.05) %>%
	pull(Gene)




pst.ord <- order(slingPseudotime(subset.sling)[,1], na.last = NA)
heatdata <- counts(sce)[, pst.ord]
heatclus <- factor(as.character(df[pst.ord,3]))


genes = rownames(heatdata)
mat = as.matrix(heatdata)
annot = as.data.frame(cbind("Int2 lineage" = identifiers[as.numeric(as.character(heatclus))+1]))
rownames(annot) = colnames(mat)

meanCtr<-function(x){
	annAll <- dimnames(x)
	means <- apply(x,1,mean,na.rm=T)
	x <- t(scale(t(x),center=means,scale=F))
	dimnames(x) <- annAll
	return(x)
}

mat.ctr = medianCtr(mat)

colors <- list('Int2 lineage' = c('24_GE_M' = "#00B9E3", '29_GE_M' = "#DB72FB", '08_GE_G1S' = "#93AA00", '33_GE_G1S' = "#FF61C3", '20_Int4' = "#00C19F", '13_Int3' = "#00BA38", '25_Int2_Zic+' = "#619CFF", '05_Int2' = "#F8766D"))
colAnn <- HeatmapAnnotation(df = annot,which = 'col',col = colors,annotation_width = unit(c(1, 4), 'cm'),gap = unit(1, 'mm'),
	annotation_legend_param = list(	
		"Int2 lineage" = list(
        		title = "Clusters",
        		at = c("24_GE_M", "29_GE_M", "08_GE_G1S", "33_GE_G1S", "20_Int4", "13_Int3", "25_Int2_Zic+", "05_Int2"),
        		labels = c("24_GE_M", "29_GE_M", "08_GE_G1S", "33_GE_G1S", "20_Int4", "13_Int3", "25_Int2_Zic+", "05_Int2")
            	)
        )
)

GOI = c("Zcchc12","Gad2","Calb2","Sst","Mcm2","Tmem130","Hap1","Ccnd1","Mki67","Ednrb","Npy","Zic1")
clustered_indices = which(genes %in% GOI)

# Plot full heatmap for supplemental figure, with genes of interest labeled
col_fun = colorRamp2(c(0, 0.75, 1.5), viridis(3))
ComplexHeatmap::Heatmap(mat.ctr,name="Expression",clustering_distance_rows = function(x) as.dist(1-cor(t(x))),clustering_method_rows="average",cluster_columns=F,column_order=rownames(slingPseudotime(subset.sling)[pst.ord,]),show_column_names = FALSE,col=col_fun,top_annotation=colAnn,right_annotation=rowAnnotation(foo = anno_mark(at = clustered_indices, labels = genes[which(genes %in% GOI)], labels_gp = gpar(fontface = "italic"))),show_row_names=F,column_names_side = "top",column_names_gp = gpar(fontsize=8),show_row_dend = FALSE,use_raster=T,raster_quality=3)

# Plot subsetted heatmap for main figure, with genes of interest plotted only
col_fun = colorRamp2(c(0, 0.5, 1), viridis(3))
ComplexHeatmap::Heatmap(mat.ctr[GOI,],name="Expression",clustering_distance_rows = function(x) as.dist(1-cor(t(x))),clustering_method_rows="average",cluster_columns=F,column_order=rownames(slingPseudotime(subset.sling)[pst.ord,]),show_column_names = FALSE,col=col_fun,top_annotation=colAnn,show_row_names=T,column_names_side = "top",column_names_gp = gpar(fontsize=8),row_names_gp = gpar(fontface = "italic"), show_row_dend = FALSE,use_raster=T,raster_quality=3)

