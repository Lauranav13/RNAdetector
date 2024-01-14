#!/usr/bin/env Rscript

#install.packages("pheatmap")
#install.packages("qvalue")
#install.packages("preprocessCore")

suppressPackageStartupMessages({
  library(optparse)
  library(rjson)
  library(dplyr)
  library(metaseqR)
  library(qvalue)
  library(ggplot2)
  library(pheatmap)
  library(preprocessCore)
})


unslash <- function(dirs) (sub("/+$", "", dirs))

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NULL, help = "config file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (opt$help) {
  print_help(opt_parser)
}

if (is.null(opt$config) || !file.exists(opt$config)) {
  print_help(opt_parser)
  stop("Config file is required!", call. = FALSE)
}

#config <- fromJSON(file = "/rnadetector/ws/storage/app/public/jobs/245/deg_config_test.json")

config <- fromJSON(file = opt$config)

if (is.null(config$description.file) || !file.exists(config$description.file)) {
  stop("Missing data description file.", call. = FALSE)
}
if (is.null(config$data.file) || !file.exists(config$data.file)) {
  stop("Missing data file.", call. = FALSE)
}
if (is.null(config$data.type)) {
  config$data.type <- "gene"
}
if (is.null(config$conditions.variables) || length(config$conditions.variables) <= 0) {
  stop("Missing condition variables.", call. = FALSE)
}
if (is.null(config$contrasts) || length(config$contrasts) <= 0) {
  stop("Missing contrasts.", call. = FALSE)
}
if (is.null(config$output.directory)) {
  stop("Missing output directory", call. = FALSE)
}
config$output.directory <- unslash(config$output.directory)

descriptions <- read.delim(config$description.file, stringsAsFactors = FALSE)
if (!("SampleName" %in% colnames(descriptions))) {
  stop("Invalid description file.", call. = FALSE)
}
data <- read.delim(config$data.file, stringsAsFactors = FALSE)
variables    <- config$conditions.variables[config$conditions.variables %in% colnames(descriptions)]

#Obtenemos nombres de condiciones, muestras y variables.
condition_names <- config$conditions.variables
sample_names <- make.names(descriptions$SampleName)
variable_names <- make.names(descriptions$SampleType)


#Construimos lista de muestras y contrastes.
samples.list <- split(sample_names, descriptions[, condition_names])

clear.names <- function (x) (gsub("\\s+", "_", x, perl = TRUE))

descriptions$condition <- clear.names(do.call(paste, 
                                              c(
                                                lapply(variables, function (v)(descriptions[[v]])), 
                                                list(sep="_"))))

contrasts.list <- sapply(config$contrasts, function(x)(paste(make.names(clear.names(x)), collapse = "_vs_")))

contrast_conditions <- make.names(config$contrasts)
contrasts <- contrasts.list

#Imprimimos por pantalla los valores de las variables para comprobar que se esta realizando bien el proceso.
cat("samples.list:\n")
for (key in names(samples.list)) {
  cat(key, ":", samples.list[[key]], "\n")
}

cat("contrasts.list:\n")
for (contrast_item in contrasts.list) {
  cat(contrast_item, "\n")
}


source.data <- data
data$gc <- 0
common.cols <- c("id", "name", "chr", "start", "end", "strand", "length", "gc")
if (all("mapped_id" %in% colnames(data))) {
  data$mapped_id <- NULL
  data <- unique(data)
}
other.cols <- setdiff(colnames(data), common.cols)
data <- data[, c(common.cols, other.cols)]

if (is.null(config$parameters)) {
  params <- list()
} else {
  params <- config$parameters
}

check.list <- function(what, defaults = NULL) {
  if (is.null(what) || !is.list(what)) {
    res <- list()
  } else {
    res <- what
  }
  if (!is.null(defaults) && is.list(defaults)) {
    res <- modifyList(defaults, res)
  }
  return (res)
}

check.vector <- function(what, default) {
  if (is.null(what) || !is.vector(what)) {
    res <- default
  } else {
    res <- what
  }
  return (res)
}

restrict.cores <- 1 / check.vector(params$num.cores, 1)
restrict.cores <- max(min(restrict.cores, 1), 1 / detectCores())
pcut <- check.vector(params$pcut, 0.05)
pcut <- max(min(pcut, 0.1), 0.01)
when.apply.filter <- check.vector(params$when.apply.filter, "prenorm")
norm.algo <- check.vector(params$norm, "edger")
norm.algo.params <- check.list(params$norm.args, get.defaults("normalization", norm.algo))
if (norm.algo == "deseq") {
  if (!is.null(norm.algo.params$locfunc) && is.character(norm.algo.params$locfunc)) {
    if (norm.algo.params$locfunc == "shorth") {
      norm.algo.params$locfunc <- genefilter::shorth
    } else {
      norm.algo.params$locfunc <- stats::median
    }
  }
}
stats.algo <- check.vector(params$stats, "limma")
stats.algo.params <- check.list(params$stats.args)
stats.algo.params <- setNames(
  lapply(stats.algo, function(a) {
    tmp <- check.list(stats.algo.params[[a]], get.defaults("statistics", a))
  }),
  stats.algo
)
default.filters <- get.defaults("gene.filter", "")
gene.filters <- check.list(params$filters)
for (f in names(gene.filters)) {
  if (is.null(gene.filters[[f]])) {
    gene.filters[[f]] <- NULL
  } else {
    gene.filters[[f]] <- modifyList(default.filters[[f]], gene.filters[[f]])
  }
}

#Si queremos realizar una normalizacion quantile-quantile, la cual no se encuentra en el propio programa, debemos descomentar los siguientes comandos que poseen #.

#data_selected <- data[, !(names(data) %in% c("id", "name", "chr", "start", "end", "strand", "length", "gc"))]
#data_normalized <- as.data.frame(normalizeQuantiles(as.matrix(data_selected)))

#data_excluded <- data[, c("id", "name", "chr", "start", "end", "strand", "length", "gc")]
#data_final <- cbind(data_excluded, data_normalized)

#Con normalización quantile-quantile, hay que sustituir data$id -> data_final$id.

suppressWarnings({
  data$id <- make.unique(data$id)
  #norm.algo.params$normalization <- "none"
  #norm.algo.params$method <- "none"
  print(norm.algo.params)
  print(stats.algo.params)
  meta <- metaseqr(
    counts = data, #_final,
    sample.list = samples.list,
    contrast = contrasts.list,
    id.col = 1,
    name.col = 2,
    gc.col = 8,
    annotation = "embedded",
    org = "custom",
    trans.level = config$data.type,
    count.type = "gene",
    when.apply.filter = when.apply.filter,
    normalization = norm.algo,
    norm.args = norm.algo.params,
    statistics = stats.algo,
    stat.args = stats.algo.params,
    pcut = pcut,
    adjust.method = check.vector(params$adjust.method, "qvalue"),
    meta.p = check.vector(params$meta.p.method, "simes"),
    gene.filters = gene.filters,
    qc.plots = c(
      "mds", "biodetection", "countsbio", "saturation", "readnoise",
      "filtered", "correl", "pairwise", "boxplot", # "lengthbias",
      "deheatmap", "volcano", # "meandiff", "meanvar", "rnacomp",
      "biodist", "venn"
    ),
    fig.format = check.vector(params$fig.formats, c("png", "pdf")),
    export.where = config$output.directory,
    export.what = c(
      "annotation", "p.value", "adj.p.value", "meta.p.value",
      "adj.meta.p.value", "fold.change", "stats"
    ),
    export.scale = c("natural", "log2"),
    export.values = c("raw", "normalized"),
    export.stats = c("mean", "median", "sd"),
    export.counts.table = TRUE,
    out.list = TRUE
  )

  save(
    source.data, meta, samples.list, contrasts.list,
    file = paste0(config$output.directory, "/data/analysis.RData")
  )
})

contrasts_to_analyze <- config$contrasts
for (contrast_pair in contrasts_to_analyze) {
   contrast_pair <- unlist(strsplit(contrast_pair, "_vs_"))
   contrast1 <- contrast_pair[1]
   contrast2 <- contrast_pair[2]
}
output_directory <- config$output.directory

######

#Obtenemos histograma de pvalue y plots de qvalue para determinar la calidad de los datos y si su correccion del pvalue es buena.

bash_command <- paste(
  "bash -c",
  shQuote(
    paste(
      "zcat",
      file.path(output_directory, "lists", paste("metaseqr_all_out_", contrast1, "_vs_", contrast2, ".txt.gz", sep = "")),
      "| cut -f 10 >",
      file.path(output_directory, "lists", paste("metaseqr_all_out_", contrast1, "_vs_", contrast2, ".txt", sep = "")),
      "&& sed -i '/^NA$/d'",
      file.path(output_directory, "lists", paste("metaseqr_all_out_", contrast1, "_vs_", contrast2, ".txt", sep = ""))
    )
  )
)
system(bash_command)

file_name <- paste(output_directory, "/lists/metaseqr_all_out_", contrast1, "_vs_", contrast2, ".txt", sep="")
p_values <- read.csv(file_name, header = TRUE)
p_values <- as.numeric(p_values$p.value_edger)

png(file = paste(output_directory, "/plots/qc/qvalue_plots_", contrast1, "_vs_", contrast2, ".png", sep=""))

qobj <- qvalue(p_values)
q_plot <- plot(qobj)
print(q_plot)
dev.off()

png(file = paste(output_directory, "/plots/qc/pvalue_hist_", contrast1, "_vs_", contrast2, ".png", sep=""))
qobj <- qvalue(p_values)
q_hist <- hist(qobj$qvalues, main = "P-value histogram", xlab = "p-values", ylab = "Frecuency")
print(q_hist)
dev.off()

#####

#MA plot PRE-NORMALIZACIÓN
#Obtenemos el MA plot para observar la distribucion de los datos prenormalizados.

bash_MA <- paste(
  "bash -c",
  shQuote(
    paste(
      "zcat",
      file.path(output_directory, "lists", paste("metaseqr_all_out_", contrast1, "_vs_", contrast2, ".txt.gz", sep = "")),
      "| cut -f 17,29 >",
      file.path(output_directory, "lists", paste("log2_fold_change_", contrast1, "_vs_", contrast2, ".txt", sep = "")),
      "&& sed -i '/^NA$/d'",
      file.path(output_directory, "lists", paste("log2_fold_change_", contrast1, "_vs_", contrast2, ".txt", sep = ""))
    )
  )
)
system(bash_MA)

file_name <- paste(output_directory, "/lists/log2_fold_change_", contrast1, "_vs_", contrast2, ".txt", sep="")
log2_change <- read.table(file_name, header = TRUE)
log2_change_1 <- as.numeric(log2_change[[paste("log2_normalized_mean_counts_", contrast1, sep = "")]])
log2_change_2 <- as.numeric(log2_change[[paste("log2_normalized_mean_counts_", contrast2, sep = "")]])

#MA PLOT CON NORMALIZACION QUANTILE-QUANTILE, otra vez deberiamos descomentar los comandos que llevan #, para modificar el codigo y que el MA plot se realizara con los datos normalizados por quantile-quantile. 
#df <- data.frame(data_normalized)

#print(colnames(df))

#counts_1 <- df[[contrast1]]
#counts_2 <- df[[contrast2]]

pdf(file = paste(output_directory, "/plots/qc/MA_prenormalization_", contrast1, "_vs_", contrast2, ".pdf", sep=""))

#A <- log((counts_1 + counts_2)/ 2)
#M <- log(counts_1) - log(counts_2)
#median_M <- median(M)

A = (log2_change_1 + log2_change_2)/2
M <- log2_change_1 - log2_change_2
median_M <- median(M)

plot(A, M, main="MA plot", xlab="A", ylab="M", col="blue", pch=20, cex=0.6)
abline(h = median_M, col = "red", lwd = 2)

legend("bottomright", legend = paste("Median:", round(median_M, 2)), col = "red", lwd = 2, bty = "n")

dev.off()

#MA plot POST-NORMALIZACIÓN
#Obtenemos el MA plot para observar la distribucion de los datos postnormalizados. Asi fue como decidimos que normalizacion se ajustaba mejor al ensayo.

bash_MA_post <- paste(
  "bash -c",
  shQuote(
    paste(
      "zcat",
      file.path(output_directory, "lists", paste("normalized_counts_table.txt.gz", sep = "")),
      "| cut -f 10,11,12 >",
      file.path(output_directory, "lists", paste("normalized_counts_table.txt", sep = "")),
      "&& sed -i '/^NA$/d'",
      file.path(output_directory, "lists", paste("normalized_counts_table.txt", sep = ""))
    )
  )
)
system(bash_MA_post)

file_name <- paste(output_directory, "/lists/normalized_counts_table.txt", sep = "")
log2_change <- read.table(file_name, header = TRUE)
log2_change_3 <- as.numeric(log2_change[[paste(contrast1)]])
log2_change_4 <- as.numeric(log2_change[[paste(contrast2)]])

pdf(file = paste(output_directory, "/plots/qc/MA_postnormalization_", contrast1, "_vs_", contrast2, ".pdf", sep=""))

A = log(log2_change_3 + log2_change_4)/2
M <- log(log2_change_3) - log(log2_change_4)

inf_indices <- which(!is.finite(M))
M[inf_indices] <- NA
median_M <- median(M, na.rm = TRUE)


plot(A, M, main="MA plot", xlab="A", ylab="M", col="blue", pch=20, cex=0.6)
abline(h = median_M, col = "red", lwd = 2)

legend("bottomright", legend = paste("Median:", round(median_M, 2)), col = "red", lwd = 2, bty = "n")

dev.off()

######

#Obtenemos la tabla descomprimida con los genes diferencialmente expresados ordenados numericamente por el log2_fold_change de mas diferencialmente expresado a menos, junto con su p_value y su FDR.


bash_log2_change <- paste(
  "bash -c",
  shQuote(
    paste(
      "zcat",
      file.path(output_directory, paste("lists/metaseqr_sig_out_", contrast1, "_vs_", contrast2, ".txt.gz", sep = "")),
      "| cut -f 1,10,11,13,17,29 >",
      file.path(output_directory, paste("lists/metaseqr_sig_out_", contrast1, "_vs_", contrast2, ".txt", sep = "")),
      "&& sed -i '/NA/d'",
      file.path(output_directory, paste("lists/metaseqr_sig_out_", contrast1, "_vs_", contrast2, ".txt", sep = ""))
    )
  )
)

system(bash_log2_change)

result <- read.table(paste(output_directory, paste("lists/metaseqr_sig_out_", contrast1, "_vs_", contrast2, ".txt", sep = ""), sep="/"), header = TRUE, sep = "\t")
result <- result[order(-abs(result[[paste("log2_normalized_fold_change_", contrast1, "_vs_", contrast2, sep = "")]])), ]

output_table <- paste(output_directory, paste("lists/metaseqr_sig_out_", contrast1, "_vs_", contrast2, ".txt", sep = ""), sep="/")
write.table(result, file = output_table, sep = "\t", row.names = FALSE)

#Eliminamos los archivos que ya no nos aportan información ni son necesarios.

file_1 <- paste(output_directory,paste("lists/log2_fold_change_", contrast1, "_vs_", contrast2, ".txt", sep = ""), sep="/")
file_2 <- paste(output_directory, paste("lists/metaseqr_all_out_", contrast1, "_vs_", contrast2, ".txt", sep = ""), sep="/")

file.remove(file_1)
file.remove(file_2)

######

#Obtenemos un heatmap con un rango de color apropiado que nos muestre los 90 genes mas diferenciados y un cluster entre ellos. 

data_matrix <- as.matrix(meta$data[[paste(contrast1,"_vs_",contrast2, sep = "")]][, c(paste("log2_normalized_mean_counts_", contrast1, sep = ""), paste("log2_normalized_mean_counts_", contrast2, sep = ""), paste("log2_normalized_fold_change_", contrast1, "_vs_", contrast2, sep = ""))])
data_matrix <- as.data.frame(data_matrix)

column_index <- paste("log2_normalized_fold_change_", contrast1, "_vs_", contrast2, sep = "")
abs_sorted_values <- data_matrix[order(-abs(data_matrix[, column_index])), ]
top_90_abs_values <- abs_sorted_values[1:min(90, nrow(abs_sorted_values)), ]

selected_data <- top_90_abs_values[, 1:2]
colnames(selected_data) <- c(contrast1, contrast2)
heatmap <- pheatmap(selected_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   show_rownames = TRUE,
                   show_colnames = TRUE, 
                   main = paste("Top 90 DEGs ", contrast1, "_vs_", contrast2, sep = ""))

file_name <- paste(output_directory, "/plots/qc/heatmap_plot_", contrast1, "_vs_", contrast2, ".png", sep = "")
png(file_name, width = 400, height = 1000, bg = "white")
print(heatmap)
dev.off()