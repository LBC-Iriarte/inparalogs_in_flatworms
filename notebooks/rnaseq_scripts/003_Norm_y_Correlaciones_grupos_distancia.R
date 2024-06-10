library(edgeR)
library(corrr)
library(ggplot2)

# setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Schistosoma")
# setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Schistosoma_Protasio2012")
setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Smansoni_SangerInstitute")

set_datos <- gsub("_all_counts.tab", "" ,list.files(path=".", pattern="*_all_counts.tab"))

# Nota, el calculo de pseudocounts falla en algunas muestras cuando se incluye los una matriz de modelado
# esto parece ser un problema de como se distribuyen los conteos en estos trabajos (no estoy eliminando gene scon 0)
# y es solo un problema con esa funciòn particular. No incluir la matriz cambia como se estiman los resultados pero el
# efecto no es del todo claro.

wanted_groups <- c("1_Close", "2_Far", "3_Diferent_Chr", "Other")
registro_correlaciones = data.frame()  
registro_uso = data.frame()  

for (data in set_datos) {
  grupos_imparalogos <- read.table("distance_corr.txt", sep = "\t")
  counts_file <- gsub(" ","", paste(data ,"_all_counts.tab"))
  group_file <- gsub(" ","", paste(data ,"_grupos.tab"))
  imparalog_file <- "Smansoni_imparalogs.ids"
  full_plot_file <- gsub(" ","", paste(data ,"correlation_by_distance_boxplot.pdf"))
  print_registro_uso <- gsub(" ","", paste(data ,"_Registro_uso.txt"))
  print_correlation_tab <- gsub(" ","", paste(data ,"_correlation_final.tab"))
  print_norm_counts <- gsub(" ","", paste(data ,"_norm_counts.tab"))  

  gene_counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names="Term_Count")

  edgeR_Input <- DGEList(counts=gene_counts)
  edgeR_Input <- calcNormFactors(edgeR_Input, method="TMM")
  edgeR_Input <- estimateCommonDisp(edgeR_Input)
  work_pseudocounts <- as.data.frame(edgeR_Input$pseudo.counts)
  
  for (c_group in wanted_groups ) {
    current_group <- gsub("3_","", gsub("2_","", gsub("1_","", c_group)))

    run_imparalog <- grupos_imparalogos[which(grupos_imparalogos$V4 == current_group),]
    Test <- grupos_imparalogos[which(grupos_imparalogos$V4 == current_group),]

    recorrer <- 1
    stop <- nrow(run_imparalog)
    
    while (recorrer <= stop ) {
      Imparalog <- run_imparalog[recorrer,1]
      gen1 <- run_imparalog[recorrer,2]
      gen2 <- run_imparalog[recorrer,3]
      distance <- run_imparalog[recorrer,5]
      
      Imparalog_pair <- gsub(" ","", paste(Imparalog,"_x_",gen1,"_x_",gen2))  
    
      counts1 <- as.numeric(work_pseudocounts[which(row.names(work_pseudocounts) %in% gen1),])
      counts2 <- as.numeric(work_pseudocounts[which(row.names(work_pseudocounts) %in% gen2),])
      
      check1=median(counts1)
      check2=median(counts2)
      
      if (check1 >= 1 && check2 >= 1) {
        use="Yes"
        pearson_res <- correlate(counts1, counts2, use = "pairwise.complete.obs", method = "pearson")
        add_information_all <- data.frame(Comparacion=Imparalog_pair, Pearson_Cor=pearson_res$x, group=c_group)
        registro_correlaciones <- rbind(as.data.frame(registro_correlaciones),as.vector(add_information_all))
      } else {
        use="No"
      }
      
      add_usage=data.frame(Imparalog=Imparalog, Gen1=gen1, Gen1=gen2, Use=use, Group=current_group, Distance=distance)
      registro_uso=rbind(as.data.frame(registro_uso),as.vector(add_usage))
      registro_uso
      recorrer <- recorrer + 1
    }
  }
  
  correlations_Close <- registro_correlaciones[which(registro_correlaciones$group == "1_Close"),]  
  correlations_Far <- registro_correlaciones[which(registro_correlaciones$group == "2_Far"),]  
  correlations_Df_Chr <- registro_correlaciones[which(registro_correlaciones$group == "3_Diferent_Chr"),]  
  
  registro_utest <- data.frame()
  current_utest <-wilcox.test(correlations_Close$Pearson_Cor, correlations_Far$Pearson_Cor)  
  add_u_test <- data.frame(Comparacion="Close_vs_Far", Method=current_utest$method , Utest_stat=current_utest$statistic, pvalue=current_utest$p.value)
  registro_utest <- rbind(as.data.frame(registro_utest),as.vector(add_u_test))
  
  current_utest <-wilcox.test(correlations_Close$Pearson_Cor, correlations_Df_Chr$Pearson_Cor)  
  add_u_test <- data.frame(Comparacion="Close_vs_Df.Chr", Method=current_utest$method , Utest_stat=current_utest$statistic, pvalue=current_utest$p.value)
  registro_utest <- rbind(as.data.frame(registro_utest),as.vector(add_u_test))
  
  current_utest <-wilcox.test(correlations_Far$Pearson_Cor, correlations_Df_Chr$Pearson_Cor)  
  add_u_test <- data.frame(Comparacion="Far_vs_Df.Chr", Method=current_utest$method , Utest_stat=current_utest$statistic, pvalue=current_utest$p.value)
  registro_utest <- rbind(as.data.frame(registro_utest),as.vector(add_u_test))
  
  registro_correlaciones <- registro_correlaciones[which(registro_correlaciones$Pearson_Cor != "NA"),]
  
  write.table(registro_utest, "Utest_Distancia.csv", sep="\t")
  write.table(registro_correlaciones, "Correlaciones_por_Distancia.csv", sep="\t")
  write.table(registro_uso, print_registro_uso, sep="\t")
  
  grafica <- ggplot(registro_correlaciones, aes(x=group, y=Pearson_Cor, fill=group )) + 
    geom_boxplot(outlier.shape = NA, show.legend = FALSE)  +
    theme_minimal() +
    labs(title="Pearson Correlation by distance", x =element_blank(), y = "Correlation") +
    scale_x_discrete(labels=c("Close", "Far", "Dif. Chr")) +
    theme(legend.title=element_blank()) +
    scale_fill_discrete(name=NULL)
  ggsave(full_plot_file,plot=grafica, width=8, height=6)
} 

