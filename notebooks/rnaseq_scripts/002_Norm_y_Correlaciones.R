library(edgeR)
library(corrr)
library(ggplot2)

setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Schistosoma")
# setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Schistosoma_Protasio2012")
# setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Smansoni_SangerInstitute")

set_datos <- gsub("_all_counts.tab", "" ,list.files(path=".", pattern="*_all_counts.tab"))

# Funcion robada para calcular la distancia de Manhatan
manhattanDistance <- function(vect1, vect2){
  dist <- abs(vect1 - vect2)
  dist <- sum(dist)
  return(dist)
}

# Nota, el calculo de pseudocounts falla en algunas muestras cuando se incluye los una matriz de modelado
# esto parece ser un problema de como se distribuyen los conteos en estos trabajos (no estoy eliminando gene scon 0)
# y es solo un problema con esa funciòn particular. No incluir la matriz cambia como se estiman los resultados pero el
# efecto no es del todo claro.

for (data in set_datos) {
  if (data == "H_microstoma" ) {
    replace_labels <- c("Egg", "Cysticercoid", "Adult")
  } else if (data == "H_microstoma_Olson2018") {
    replace_labels <- c("Neck-Scolex", "Mid-Strobila", "Posterior-Strobila")
  } else if (data == "Schistosoma") {
    replace_labels <- c("D6", "D13", "D17", "D21", "D28", "D35")
  } else if (data == "Schistosoma_Protasio2012") {
    replace_labels <- c("Cercariae", "Schistosomulae 3hr", "Schistosomulae 24hr")
  } else if (data == "Smansoni_SangerInstitute") {
    replace_labels <- c("Miracidia", "First stage sporocyst", "Second stage sporocyst")    
  }
  
  counts_file <- gsub(" ","", paste(data ,"_all_counts.tab"))
  group_file <- gsub(" ","", paste(data ,"_grupos.tab"))
  imparalog_file <- "Smansoni_imparalogs.ids"
  
  print_correlation_tab <- gsub(" ","", paste(data ,"_correlation_final.tab"))
  print_norm_counts <- gsub(" ","", paste(data ,"_norm_counts.tab"))  

  gene_counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names="Term_Count")
  group_list <- read.table(group_file, header = TRUE, sep = "\t")
  imparalog_list <- read.table(imparalog_file, header = FALSE, sep = "\t", row.names="V1")
  
  edgeR_Input <- DGEList(counts=gene_counts)
  edgeR_Input <- calcNormFactors(edgeR_Input, method="TMM")
  edgeR_Input <- estimateCommonDisp(edgeR_Input)
  work_pseudocounts <- as.data.frame(edgeR_Input$pseudo.counts)
  
  write.table(work_pseudocounts,print_norm_counts, sep="\t")
  
  run_imparalog <- row.names(imparalog_list)
  run_group <- unique(group_list$Grupo)
  
  correlaciones_resultados_finales = data.frame()
  correlaciones_chisquare_finales = data.frame()
  
  for (imp in run_imparalog) {
    print(imp)
    imp_title=paste("Monophyletic group:", imp)

    save_counts_file <- gsub(" ","", paste(imp ,"_imparalog_pseudocounts.txt"))
    median_plot_file <- gsub(" ","", paste(imp ,"_median_plot.pdf"))
    full_plot_file <- gsub(" ","", paste(imp ,"_full_plot.pdf"))
    
    string1 <- gsub(",", "", imparalog_list[imp,])
    extract_genes <- strsplit(string1, " ")[[1]]
    temp_data <- work_pseudocounts[which(row.names(work_pseudocounts) %in% extract_genes),]
    
    current_imparalog_median_data = data.frame()
    current_imparalog_all_data = data.frame()
    for (gru in run_group) {
      group_genes <- group_list[group_list$Grupo == gru, 1 ]
      tem2 <- temp_data[,which(names(temp_data) %in% group_genes)]
      
      for (gen in extract_genes) {
        # Calculate medians
        group_median <- median(as.numeric(as.vector(tem2[which(row.names(tem2) %in% gen),])))
        add_information_median <- data.frame(Gene_id=gen, Median=group_median, group=gru)
        current_imparalog_median_data <- rbind(as.data.frame(current_imparalog_median_data),as.vector(add_information_median))
        as.vector(tem2[which(row.names(tem2) %in% gen),])
        # Hacer grafica on todos los conteos todos
        all_gene_counts <- as.numeric(as.vector(tem2[which(row.names(tem2) %in% gen),]))
        rep_elementos_conteo <- length(all_gene_counts)
        print_gene <- rep.int(gen, rep_elementos_conteo)
        print_gru <- rep.int(gru, rep_elementos_conteo)
        
        add_information_all <- data.frame(Gene_id=rep.int(gen, rep_elementos_conteo), Pseudo_counts=all_gene_counts, group=rep.int(gru, rep_elementos_conteo))
        current_imparalog_all_data <- rbind(as.data.frame(current_imparalog_all_data),as.vector(add_information_all))
      }
    }
    
    grafica <- ggplot(current_imparalog_median_data, aes(x=group, y=Median, color=Gene_id, group =Gene_id)) + 
      geom_point(data = current_imparalog_all_data, aes(x=group, y=Pseudo_counts)) +
      geom_line() +
      theme_minimal() +
      labs(title=imp_title, x ="Group", y = "TMM") +
      geom_hline(yintercept=1, linetype="dashed", color = "red") +
      theme(legend.title=element_blank()) +
      scale_x_discrete(labels=replace_labels) +
      scale_fill_discrete(name=NULL)
    ggsave(full_plot_file,plot=grafica, width=8, height=6)

    
    ggplot(current_imparalog_median_data, aes(x=group, y=Median, color=Gene_id, group =Gene_id)) + 
      geom_point(data = current_imparalog_all_data, aes(x=group, y=Pseudo_counts)) +
      geom_line() +
      theme_minimal() +
      labs(title=imp_title, x ="Group", y = "TMM") +
      geom_hline(yintercept=1, linetype="dashed", color = "red") +
      theme(legend.title=element_blank()) +
      scale_x_discrete(labels=replace_labels) +
      scale_fill_discrete(name=NULL)
    
    # Calculate all correlations
    correlaciones_completadas <- as.data.frame("DoNothing")
    for (gen1 in extract_genes) {
      for (gen2 in extract_genes) {
        if (gen1 != gen2) {
          temp_registro <- paste(sort(c(gen1, gen2)),collapse=" ")
          if ( ! temp_registro %in% as.character(as.matrix(correlaciones_completadas)) ) {
            correlaciones_completadas <- rbind(as.data.frame(correlaciones_completadas), as.vector(temp_registro))     
            data_gen1 <- current_imparalog_median_data[current_imparalog_median_data$Gene_id == gen1,2]
            data_gen2 <- current_imparalog_median_data[current_imparalog_median_data$Gene_id == gen2,2]
            
            # Calcular correlaciones
            pearson_res <- correlate(data_gen1, data_gen2, use = "pairwise.complete.obs", method = "pearson")
            spearman_res <- correlate(data_gen1, data_gen2, use = "pairwise.complete.obs", method = "spearman")
            
            print_manhattanDistance_res <- manhattanDistance(data_gen1, data_gen2)

            print_resultados <- data.frame(Inparalogs_ID=imp, GenA=gen1, GenB=gen2, Pearson=pearson_res$x, Speerman=spearman_res$x, ManhattanDistance = print_manhattanDistance_res)
            correlaciones_resultados_finales <- rbind(as.data.frame(correlaciones_resultados_finales), as.vector(print_resultados))

          }
        }
      }
    }
  }
  write.table(correlaciones_resultados_finales,print_correlation_tab, sep="\t")
}



