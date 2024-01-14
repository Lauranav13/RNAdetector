library(pheatmap)
library(gplots)

data_matrix <- read.csv("/home/vant/pLs20_Experiments/DEG_Analysis/ARCHIVOS A ENTREGAR TFM/Imágenes/Resultados/PKS11_vs_PKS14/pks11_vs_pks14.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

print(colnames(data_matrix))

data_matrix[, c("log2_normalized_fold_change_pks11_vs_pks14", "log2_normalized_pks11", "log2_normalized_pks14")] <- lapply(
  data_matrix[, c("log2_normalized_fold_change_pks11_vs_pks14", "log2_normalized_pks11", "log2_normalized_pks14")],
  function(x) as.numeric(gsub(",", ".", x))
)

column_index <- paste("log2_normalized_fold_change_pks11_vs_pks14", sep = "")
abs_sorted_values <- data_matrix[order(-abs(data_matrix[, column_index])), ]
top_90_abs_values <- abs_sorted_values[1:min(90, nrow(abs_sorted_values)), ]

selected_data <- top_90_abs_values[, c(
  paste("log2_normalized_pks11", sep = ""),
  paste("log2_normalized_pks14", sep = "")
)]
colnames(selected_data) <- c("pks11", "pks14")
colnames(selected_data)

selected_data <- selected_data[, c("pks11", "pks14")]

heatmap <- pheatmap(selected_data, 
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE, 
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    labels_row = top_90_abs_values$GeneID,  
                    main = paste("Top 90 DEGs pks11_vs_pks14 without pLS20", sep = ""))

file_name <- "/home/vant/pLs20_Experiments/DEG_Analysis/ARCHIVOS A ENTREGAR TFM/Imágenes/Resultados/PKS11_vs_PKS14/heatmap_withoutpLS20_Q-Q.png"
png(file_name, width = 400, height = 1000, bg = "white")
print(heatmap)
dev.off()


