library(xml2)
library(jsonlite)
library(dplyr)
library(purrr)

#Introducimos el archivo GTF utilizado en el análisis de expresión diferencial de genes para sacar los ID con su respectivo GeneID.
data <- readLines("/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168.gtf")

gene_ids <- character()
gene_names <- character()

#Iteramos sobre las líneas del archivo.
for (line in data) {
  #Buscamos GENEID:
  gene_id_match <- regmatches(line, regexpr('GeneID:([0-9]+)', line))
  
  if (length(gene_id_match) > 0) {
    gene_ids <- c(gene_ids, gene_id_match)
  }
  
  #Buscamos sobre las líneas gene:
  gene_name_match <- regmatches(line, regexpr('gene "([^"]+)"', line))
  
  if (length(gene_name_match) > 0) {
    gene_names <- c(gene_names, gene_name_match)
  }
}

#Creamos un dataframe con los resultados.
result_df <- data.frame(ID = gene_ids, GeneID = gene_names)

#Guardamos el resultado en un nuevo archivo GTF.
write.table(result_df, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168_1.gtf", sep = "\t", row.names = FALSE, quote = FALSE)

data <- readLines("/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168_1.gtf")

#Depuramos el vector de caracteres.
data <- gsub("GeneID:", "", data)
data <- gsub('gene "', '', data)
data <- gsub('"', '', data)

#Escribimos el vector de caracteres depurado en un nuevo archivo.
writeLines(data, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168_cleaned.gtf")

result_df <- read.table("/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168_cleaned.gtf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Ordenamos el dataframe por la columna "GeneID".
result_df <- result_df[order(result_df$GeneID), ]

#Guardamos el dataframe ordenado en un nuevo archivo.
write.table(result_df, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168_sorted.gtf", sep = "\t", row.names = FALSE, quote = FALSE)

#Metemos el resultado en una variable para poder operar con ella.
geneid_data <- read.table("/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/Bsub168_sorted.gtf", header = TRUE, stringsAsFactors = FALSE)
#Después cogemos esos ID para meterlos en otro archivo y usarlo posteriormente en la hipergeométrica.
write.table(geneid_data$ID, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/all_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Ahora tenemos todos los genes de B. subtilis 168 con su GeneID y su ID. El siguiente paso es filtrar los genes por los que seleccionamos anteriormente mediante el DEG.

#Este tiene que ser el archivo resultante del DEG, con el GeneID, el Qvalue y el Log2FoldChange. Esto se saca seleccionando esas columnas en el metaseqr_sig_out.txt y transformandolo en un txt.
genes_to_filter <- read.table("/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/genes.txt", header = FALSE, col.names = c("GeneID", "Log2FoldChange"))
genes_to_filter <- genes_to_filter[order(genes_to_filter$GeneID), ]
write.table(genes_to_filter, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Se une con el archivo creado para filtrar los genes resultantes del DEG y obetener un archivo con el GeneID, el ID y el Log2FoldChange.
filtered_data <- inner_join(geneid_data, genes_to_filter, by = "GeneID")
str(filtered_data)

#Por último se guarda en una tabla.
write.table(filtered_data, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/filtered_geneid.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Después cogemos esos ID diferencialmente expresados para meterlos en otro archivo y usarlo posteriormente en la hipergeométrica.
write.table(filtered_data$ID, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/differentially_expressed_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#El JSON esta sacado del NCBI, es una base de datos de rutas metabólicas de Bacillus subtilis.
json_data <- fromJSON("/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/PubChem_pathway_text_bacillus_subtilis.json")

#Cogemos los ID de la base de datos.
json_data$geneid <- lapply(strsplit(as.character(json_data$geneids), " "), as.integer)
geneids <- json_data$geneids

data <- data.frame(ID = I(geneids), stringsAsFactors = FALSE)

#Cogemos los pathways de la base de datos.
pathways <- json_data$pwacc

data_pathway <- data.frame(PATHWAY = pathways, stringsAsFactors = FALSE)

if (nrow(data) == nrow(data_pathway)) {
  merged_data <- cbind(data, data_pathway)
  print(merged_data)
}
str(merged_data)

write.table(merged_data, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/ID_PATHWAY.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Guardamos una tabla con los ID y los PATHWAY en la que los ID no estén escritos en forma de vector y sin los NULL.
merged_table <- merged_data
merged_table$ID <- gsub("[\"c()]", "", merged_table$ID)
merged_table <- merged_table[!sapply(merged_table$ID, function(x) is.null(x) || x == ""), ]
merged_table <- merged_table[merged_table$ID != "NULL", ]

write.table(merged_table, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/ID_PATHWAY_TABLE.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Después cogemos esos ID para meterlos en otro archivo y usarlo posteriormente para ver los IDs coincidentes en bash.
write.table(merged_table$ID, "/home/vant/pLs20_Experiments/DEG_Analysis/Enrichment_analysis/ID_TABLE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
