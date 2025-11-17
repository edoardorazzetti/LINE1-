library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
#a <- read.delim(“0_Input_DEG_pt_coding.tsv”, row.names=1)
colorh<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(256))
#class.c<-colorRampPalette(brewer.pal(12,“Paired”))(12)
colclass <- c("Adjacent"="orange2", "Tumor"="red3")
isoforms_scale <- scale(log10(t(isoforms+1)))
column_ha = HeatmapAnnotation(Class=metadata$Tissue,
                                             col = list(Class = colclass),
                                             border = TRUE,
                                             simple_anno_size = unit(0.3, "cm"),
                                             annotation_name_gp = gpar(fontsize =9),
                                             annotation_legend_param = list(border="black"))
Heatmap(t(isoforms_scale),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 2),
  cluster_columns = T,
  cluster_rows = T,
  name = "Z-score",
  border = "black",
  column_names_rot = 45,
  col=colorRamp2(breaks=c(-2,-1,0,1,2), colors=c(colorh[1], colorh[64], "white",colorh[192], colorh[255])),
                       heatmap_legend_param = list(border = "black"), width = 50,
                       #clustering_method_rows = “ward.D2”,
                       #clustering_method_columns = “ward.D2",
        top_annotation = column_ha)
#112 pazienti totali, considerando anche i polipi. 90 casi di crc, 22 di polipi








