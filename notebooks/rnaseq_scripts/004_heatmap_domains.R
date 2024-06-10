setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Heatmap_Files/Una_especie")
setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Heatmap_Files/Dos_Especies")

library(pheatmap)   
library(RColorBrewer)

set_datos <- gsub("_work_norm.in", "" ,list.files(path=".", pattern="*_work_norm.in"))
# breaksList = seq(0.0, 1, by = 0.1) 

for (data in set_datos) {
  data_file <- gsub(" ","", paste(data ,"_work_norm.in"))
  heatmap_file <- gsub(" ","", paste(data ,"_heatmap.pdf"))
  print_domain <- gsub(" ","", paste(data ,"_domain_list.txt"))
  
  title <- paste("Interpro Signatures for GOterm:" , gsub("_", ":",gsub("_heatmap_.*", "" ,data)))
  data_for_heatmap <- read.table(data_file, header = TRUE, row.names = "Species")
  data_for_heatmap <- replace(data_for_heatmap, data_for_heatmap>0, 1) 
  
  domains <- as.data.frame(colnames(data_for_heatmap))
  colnames(domains) <- c("Domains") 

  # data_for_heatmap <- replace(data_for_heatmap, data_for_heatmap==0, NA) 
    
  write.table(domains, print_domain)
  pdf(file=heatmap_file, width=80, height=20 )
                                                                                             
    pheatmap(data_for_heatmap,main = title, fontsize=10, cellwidth=12, cellheight=12, legend = FALSE , cluster_rows= FALSE, color=colorRampPalette(c("blue3", "red3"))(2))
    # pheatmap(data_for_heatmap,fontsize=10, cellwidth=12, cellheight=12, cluster_rows= FALSE, breaks = breaksList, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
    # color=colorRampPalette(c("ivory", "blue3"))(2) Esto reduce la paleta a 2 colores
  dev.off()
}


setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos")
library(pheatmap)   
library(RColorBrewer)
data_for_heatmap <- read.table("MeropsFamiliesXSp.tab", header = TRUE, row.names = "Orga")  
data_for_heatmap <-as.data.frame(t(apply(data_for_heatmap,1, function(x) x/sum(x))))

pdf(file="Test_Protease_heatmap_V6.pdf", width=25, height=12 )
  pheatmap(data_for_heatmap,fontsize=10, cellheight=12, cluster_rows= FALSE, cluster_cols = FALSE, na_col="white", border_color = NA, color=colorRampPalette(c("white", "red4"))(10))
dev.off()


data_for_heatmap <- replace(data_for_heatmap, data_for_heatmap==0, NA) 
data_for_heatmap <- log10(data_for_heatmap)

pdf(file="Test_Protease_heatmap_Log10_V4.pdf", width=25, height=12 )
  pheatmap(data_for_heatmap,fontsize=10, cellheight=12, cluster_rows= FALSE, cluster_cols = FALSE, na_col="white", border_color = NA, color=colorRampPalette(c("white", "red4"))(10))
dev.off()

