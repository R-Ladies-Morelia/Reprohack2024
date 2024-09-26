# Título: Análisis de genomas bacterianos evolucionados en ausencia de fagos
# Autora: Reena Debray
# Fecha: 1 de febrero de 2022
# Modificado: Amairani Cancino Bello (04/Sep/2024)

### Funciones
## Esta función utiliza la salida del programa breseq para construir una matriz de genes que fueron mutados en cada población durante la evolución experimental.
### Toma datos en la siguiente forma: "breseq", un data frame con una entrada para cada mutación en cada población en cada momento y su frecuencia alélica correspondiente.
### "gene_list", una lista de todos los genes observados en la población o subpoblación de interés.
### "sample_list", una lista de todas las muestras en la población o subpoblación de interés.
### AF_min, una frecuencia alélica mínima para su consideración. He establecido AF_min=0 (sin mínimo), aunque hay que notar que breseq solo devuelve polimorfismos con una frecuencia poblacional de 0.05 o mayor.

# --- paquetes ----
library(readxl)
library(ggplot2)
install.packages("viridis")
library(viridis)

#--- Especificar la ruta completa del directorio ----
breseq <- "C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables"

#--- Generar matriz de genes anotada ---- 
#' @title Generar matriz de genes anotada
#' @description Crea una matriz anotada con genes y muestras utilizando los datos de mutación de breseq.
#' @param breseq Data frame que contiene los datos de mutación de breseq.
#' @param gene_list Lista de genes a incluir en la matriz.
#' @param sample_list Lista de muestras a incluir en la matriz.
#' @param AF_min Frecuencia alélica mínima para filtrado.
#' @return Una matriz anotada con el número de sitios de mutación por gen por muestra.
#' @export
gene_matrix_annotated <- function(breseq, gene_list, sample_list, AF_min) {
  LOR_gene_matrix <- matrix(ncol = length(gene_list), nrow = length(sample_list))
  colnames(LOR_gene_matrix) = gene_list
  rownames(LOR_gene_matrix) = sample_list
  
  for (gene in gene_list) {
    for (sample in sample_list) {
      if (gene %in% breseq[breseq$Sample == sample & breseq$freq >= AF_min, "gene_2"]) {
        num_sites <- length(unique(breseq[breseq$Sample == sample & breseq$gene_2 == gene & !is.na(breseq$gene_2), "position"]))
        LOR_gene_matrix[sample, gene] = num_sites
      } else {
        LOR_gene_matrix[sample, gene] = 0
      }
    }
  }
  return(LOR_gene_matrix)
}

#--- Leer y filtrar la salida de breseq ----
breseq_annotated <- data.frame()
filenames = list.files("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables")
for (file in filenames) {
  ann <- data.frame(read_excel(paste("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables", file, sep = "/")))
  ann$Line <- unlist(strsplit(file, split = "_"))[1]
  ann$Passage <- unlist(strsplit(file, split = "_"))[2]
  ann$Sample <- paste(ann$Line, ann$Passage, sep = "_")
  if (ncol(ann) > 10) {
    ann <- ann[, -c(8:12)] # algunos archivos tienen una columna extra de notas para anotaciones sin evidencia (*) resultando en espacios en las columnas 8-12
  }
  breseq_annotated <- rbind(breseq_annotated, ann)
}
dim(breseq_annotated)

#--- Eliminar evidencia no asignada ----
breseq_annotated <- breseq_annotated[!breseq_annotated$evidence %in% c(NA,"Unassigned missing coverage evidence","*","Unassigned new junction evidence","?"),]
dim(breseq_annotated)

# --- Eliminar sitios que difieren del referencia ----
## Eliminar cualquier sitio que difiera del referencia en ANCDC3000
ANCDC3000_sites <- breseq_annotated[breseq_annotated$Line == "ANCDC3000", "position"]
breseq_annotated <- breseq_annotated[!breseq_annotated$position %in% ANCDC3000_sites,]
dim(breseq_annotated)

# --- Eliminar sitios fijos ----
## Eliminar sitios que están fijos en cada línea
tab <- data.frame(table(breseq_annotated[breseq_annotated$freq == 1, "Line"], breseq_annotated[breseq_annotated$freq == 1, "position"])) # tabla de sitios fijos
tab <- tab[tab$Freq > 0,] # no considerar combinaciones de línea-mutación que nunca ocurrieron
mut_tab <- data.frame(table(tab$Var2)) # totalizar el número de líneas en las que la mutación ocurrió alguna vez
sites_to_remove <- as.character(mut_tab[mut_tab$Freq == length(unique(breseq_annotated$Line)), "Var1"])
breseq_annotated <- breseq_annotated[!breseq_annotated$position %in% sites_to_remove,]
dim(breseq_annotated)

# --- Restructurar las mutaciones de desplazamiento de marco ----
## Restructurar las mutaciones de desplazamiento de marco para ser consistentes con las mutaciones puntuales
for (i in seq(1, nrow(breseq_annotated))) {
  pos <- breseq_annotated[i, "position"]
  # encontrar mutaciones con un ":"
  if (grepl(":", pos, fixed = TRUE)) {
    pos <- unlist(strsplit(pos, split = ":"))[1]
    pos <- paste(unlist(strsplit(pos, split = ",")), collapse = "")
    breseq_annotated[i, "position"] <- pos
  }
}
breseq_annotated$position <- as.numeric(breseq_annotated$position)

# --- Restructurar columna de genes ----
## Restructurar la columna de genes para eliminar la dirección
for (i in seq(1, nrow(breseq_annotated))) {
  breseq_annotated[i, "gene_2"] <- substr(breseq_annotated[i, "gene"], 1, nchar(breseq_annotated[i, "gene"]) - 2)
}

# --- Filtrar mutaciones con ventana deslizante ----
# Escanear una ventana deslizante de 50 pb dentro de cada línea/punto en el tiempo. Si una mutación tiene más de N vecinos (incluyendo llamadas repetidas en la misma posición), eliminarla.
breseq_annotated_filtered <- data.frame()
N = 4

for (sample in unique(breseq_annotated$Sample)) {
  pre_filter_POS <- sort(as.numeric(breseq_annotated[breseq_annotated$Sample == sample, "position"]))
  sites_to_remove <- c()
  for (pos in pre_filter_POS) {
    # ventana deslizante
    min = pos - 50
    while (min < pos) {
      max <- min + 50
      if (length(pre_filter_POS[pre_filter_POS >= min & pre_filter_POS < max]) > N) {
        sites_to_remove <- unique(c(sites_to_remove, pre_filter_POS[pre_filter_POS >= min & pre_filter_POS < max])) # si una ventana contiene más de 3 mutaciones, eliminar todas ellas
      }
      min <- min + 1
    }
  }
  
  ### filtrar sitios
  post_filter <- breseq_annotated[breseq_annotated$Sample == sample & !breseq_annotated$position %in% sites_to_remove,]
  breseq_annotated_filtered <- rbind(breseq_annotated_filtered, post_filter)
}

# Imprimir las dimensiones de los datos filtrados
dim(breseq_annotated_filtered)

## Resumir anotaciones como NS, S, pseudogen o intergénico
for (i in seq(1, nrow(breseq_annotated_filtered))) {
  anno <- substr(breseq_annotated_filtered[i, "annotation"], 1, 5)
  if (anno == "pseud") {
    breseq_annotated_filtered[i, "Type"] <- "Pseudogene"
  } else if (anno == "inter") {
    breseq_annotated_filtered[i, "Type"] <- "Intergenic"
  } else {
    AA1 <- substr(anno, 1, 1)
    AA2 <- substr(anno, 5, 5)
    if (AA1 == AA2) {
      breseq_annotated_filtered[i, "Type"] <- "S"
    } else {
      breseq_annotated_filtered[i, "Type"] <- "NS"
    }
  }
}

## Contar el número de sitios fijos y polimorfismos en cada población en el paso 12
### Mutaciones fijas (todos los tipos)
### Excluir mutaciones de resistencia del conteo (porque se adquirieron antes de la evolución experimental)
Costs_of_res <- read_excel("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Costs_of_res.xlsx")
res_mutations <- unique(Costs_of_res$Position)
a <- aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage == "P12" & breseq_annotated_filtered$Type %in% c("NS", "S") & !breseq_annotated_filtered$position %in% c(res_mutations), "freq"], 
               by = list(breseq_annotated_filtered[breseq_annotated_filtered$Passage == "P12" & breseq_annotated_filtered$Type %in% c("NS", "S") & !breseq_annotated_filtered$position %in% c(res_mutations), "Line"]), 
               FUN = function(x) { length(x[x == 1]) })
mean(a$x)
sd(a$x)
t.test(x ~ (substr(a$Group.1, 1, 3) == "ANC"), a, var.equal = T)

### Polimorfismos (todos los tipos)
p <- aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage == "P12" & breseq_annotated_filtered$Type %in% c("NS", "S") & !breseq_annotated_filtered$position %in% c(res_mutations), "freq"], 
               by = list(breseq_annotated_filtered[breseq_annotated_filtered$Passage == "P12" & breseq_annotated_filtered$Type %in% c("NS", "S") & !breseq_annotated_filtered$position %in% c(res_mutations), "Line"]), 
               FUN = function(x) { length(x[x < 1]) })
mean(p$x)
sd(p$x) 
t.test(x ~ (substr(p$Group.1, 1, 3) == "ANC"), p, var.equal = T)

### Mutaciones fijas en poblaciones mutadoras (dentro de genes)
mut <- c("FMS13", "MS1", "MS10", "MS15", "QAC4")
p[p$Group.1 %in% mut, "Mutator"] <- "Y"
p[!p$Group.1 %in% mut, "Mutator"] <- "N"
t.test(x ~ Mutator, p, alternative = "less", var.equal = T)

#--- Figura 4: Paralelismo genómico a través de poblaciones----
### Generar nombres de poblaciones iniciadas con bacterias resistentes a fagos solamente (excluir controles ancestrales)
n <- unique(breseq_annotated_filtered[substr(breseq_annotated_filtered$Line, 1, 3) != "ANC", "Line"])

### Formar un data frame de mutaciones no sinónimas no fijas en Psg 0 solamente
breseq_NS_noP0 <- data.frame(matrix(nrow = 0, ncol = ncol(breseq_annotated_filtered)))
for (line in n) {
  P0_mutations <- breseq_annotated_filtered[breseq_annotated_filtered$Line == line & breseq_annotated_filtered$Passage == "P0" & breseq_annotated_filtered$freq == 1, "position"]
  ### filas correspondientes a mutaciones P2/P12 no fijas en P0 de esta muestra
  sample_data <- breseq_annotated_filtered[breseq_annotated_filtered$Line == line & breseq_annotated_filtered$Passage != "P0" & !breseq_annotated_filtered$position %in% P0_mutations & breseq_annotated_filtered$Type == "NS", ]
  
  breseq_NS_noP0 <- rbind(breseq_NS_noP0, sample_data)
}

### Hacer una matriz de genes utilizando solo mutaciones no fijas en P0
gene_list <- sort(unique(breseq_NS_noP0$"gene_2"))
sample_list <- unique(breseq_NS_noP0$Sample)
LOR_gene_matrix_NS_noP0 <- gene_matrix_annotated(breseq_NS_noP0, gene_list, sample_list, 0)

### Calcular coeficientes de similitud de Sørensen-Dice para cada par de muestras
all_pairs <- combn(sample_list, 2)
sim_coef_NS_noP0 <- data.frame(matrix(nrow = 0, ncol = 3))

for (i in seq(1, ncol(all_pairs))) {
  sample1 <- all_pairs[1, i]
  sample2 <- all_pairs[2, i]
  
  ### Nombres de genes con cualquier mutación del genotipo ancestral (cualquier número de mutaciones independientes por gen)
  sample1_genes <- names(LOR_gene_matrix_NS_noP0[sample1, apply(data.frame(LOR_gene_matrix_NS_noP0[sample1, ]), 1, function(x) { x > 0 })])
  sample2_genes <- names(LOR_gene_matrix_NS_noP0[sample2, apply(data.frame(LOR_gene_matrix_NS_noP0[sample2, ]), 1, function(x) { x > 0 })])
  
  ### Calcular el grado de superposición
  sim_coef <- (2 * length(intersect(sample1_genes, sample2_genes))) / (length(sample1_genes) + length(sample2_genes))
  sim_coef_NS_noP0 <- rbind(sim_coef_NS_noP0, c(sample1, sample2, sim_coef))
}
colnames(sim_coef_NS_noP0) <- c("sample1", "sample2", "sim_coef")
sim_coef_NS_noP0$sim_coef <- as.numeric(sim_coef_NS_noP0$sim_coef)

### Anotar los genes de resistencia de las poblaciones
for (i in seq(1, nrow(sim_coef_NS_noP0))) {
  sim_coef_NS_noP0[i, "line1"] <- unlist(strsplit(sim_coef_NS_noP0[i, "sample1"], split = "_"))[1]
  sim_coef_NS_noP0[i, "passage1"] <- unlist(strsplit(sim_coef_NS_noP0[i, "sample1"], split = "_"))[2]
  sim_coef_NS_noP0[i, "line2"] <- unlist(strsplit(sim_coef_NS_noP0[i, "sample2"], split = "_"))[1]
  sim_coef_NS_noP0[i, "passage2"] <- unlist(strsplit(sim_coef_NS_noP0[i, "sample2"], split = "_"))[2]
  
  sim_coef_NS_noP0[i, "gene1"] <- as.character(unique(Costs_of_res[Costs_of_res$Population == sim_coef_NS_noP0[i, "line1"], "Gene"]))
  sim_coef_NS_noP0[i, "gene2"] <- as.character(unique(Costs_of_res[Costs_of_res$Population == sim_coef_NS_noP0[i, "line2"], "Gene"]))
}

# Anotar si los pares tienen los mismos o diferentes genes de resistencia
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1 == sim_coef_NS_noP0$gene2, "same_gene"] <- "same"
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1 != sim_coef_NS_noP0$gene2, "same_gene"] <- "different"

##--- Prueba de aleatorización (Día 6)----
### Construir estadístico de prueba
tmp_day6 <- sim_coef_NS_noP0[sim_coef_NS_noP0$passage1 == "P2" & sim_coef_NS_noP0$passage2 == "P2" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2), ]
test_same_P2_day6 <- tmp_day6[tmp_day6$same_gene == "same", "sim_coef"]
test_diff_P2_day6 <- tmp_day6[tmp_day6$same_gene == "different", "sim_coef"]
test_gene_P2_day6 <- mean(test_same_P2_day6) - mean(test_diff_P2_day6)

### Aleatorización y re-muestreo
perm_gene_P2_day6 <- c()
set.seed(123)
for (i in seq(1, 10000)) {
  tmp_day6$same_gene <- sample(tmp_day6$same_gene, replace = F)
  perm_same_P2 <- mean(tmp_day6[tmp_day6$same_gene == "same", "sim_coef"])
  perm_diff_P2 <- mean(tmp_day6[tmp_day6$same_gene == "different", "sim_coef"])
  perm_gene_P2_day6 <- c(perm_gene_P2_day6, perm_same_P2 - perm_diff_P2)
}

### p-valor (proporción de permutaciones más extremas que el valor observado)
print(length(perm_gene_P2_day6[perm_gene_P2_day6 >= test_gene_P2_day6]) / length(perm_gene_P2_day6))

#--- Modificado para agregar el valor p----

# Preparar datos
test_gene_df_day6 <- data.frame(
  sim_coef = c(test_same_P2_day6, test_diff_P2_day6),
  label = c(rep("Mismo\ngen de\nresistencia", length(test_same_P2_day6)), rep("Diferente\ngen de\nresistencia", length(test_diff_P2_day6)))
)
colnames(test_gene_df_day6) <- c("sim_coef", "label")
test_gene_df_day6$label <- factor(test_gene_df_day6$label, levels = c("Mismo\ngen de\nresistencia", "Diferente\ngen de\nresistencia"))

# Calcular el valor p del día 6
p_value_P2_day6 <- length(perm_gene_P2_day6[perm_gene_P2_day6 >= test_gene_P2_day6]) / length(perm_gene_P2_day6)
print(p_value_P2_day6)

# Configurar anotación del valor p
if (p_value_P2_day6 < 0.001) {
  annotation_day6 <- "p < 0.001***"
} else {
  annotation_day6 <- paste0("p = ", format(p_value_P2_day6, digits = 3))
}

#--- Fig. 4A ---- 
# Creación y visualización del gráfico
p_day6 <- ggplot(test_gene_df_day6, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Superposición en mutaciones adquiridas\n(Similitud de Sørensen-Dice por pares)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55) +
  annotate("text", x = 1.5, y = 0.55, label = annotation_day6, size = 6, hjust = 0.5, vjust = 1.5, color = "black")

print(p_day6)

#--- Fig. 4B ---- 
## Identificar qué genes fueron mutados en qué poblaciones
top_genes_in_study_day6 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "Line"]))
colnames(top_genes_in_study_day6) <- c("description", "Line", "hits")
top_genes_day6 <- names(tail(sort(table(top_genes_in_study_day6[top_genes_in_study_day6$hits > 0, "description"])), 20))

### Construir un data frame de mutaciones que aparecieron para el día 6, excluyendo mutaciones de resistencia a fagos o cualquier otra mutación fijada antes del inicio del experimento
breseq_gene_analysis_day6 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "Line"]))
colnames(breseq_gene_analysis_day6) <- c("gene", "Line", "hits")

# Definir los vectores o listas de valores correspondientes a cada gen de resistencia
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  

### Anotar cada población por su gen de resistencia
breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% glycoside, "res_gene"] <- "glycoside"

# Eliminar filas donde res_gene es NA
breseq_gene_analysis_day6 <- breseq_gene_analysis_day6[!is.na(breseq_gene_analysis_day6$res_gene), ]

### Calcular el porcentaje de poblaciones con cada gen de resistencia que mutó en cada uno de los otros genes
breseq_gene_props_day6 <- data.frame(matrix(nrow = 0, ncol = 3))
l <- list(rfbA, glycosyl, glycoside)
names(l) <- c("mutantes de rfbA", "mutantes de PSPTO_4988", "mutantes de PSPTO_4991")

# Bucle a través de genes y calcular proporciones de mutación
for (gene in top_genes_day6) {
  for (j in seq_along(l)) {
    res_gene <- l[[j]]
    prop <- nrow(breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% res_gene & 
                                             breseq_gene_analysis_day6$gene == gene & 
                                             breseq_gene_analysis_day6$hits > 0, ]) / 
      nrow(breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% res_gene & 
                                       breseq_gene_analysis_day6$gene == gene, ])
    breseq_gene_props_day6 <- rbind(breseq_gene_props_day6, c(gene, names(l)[j], prop))
  }
}

# Ordenar la columna 'res_gene' para que "mutantes de rfbA" aparezca primero, seguido de "mutantes de PSPTO_4988" y luego "mutantes de PSPTO_4991"
colnames(breseq_gene_props_day6)
colnames(breseq_gene_props_day6) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day6$res_gene <- factor(breseq_gene_props_day6$res_gene, levels = c("mutantes de PSPTO_4991", "mutantes de PSPTO_4988", "mutantes de rfbA"))

# Asignar nombres de columnas
colnames(breseq_gene_props_day6) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day6$percent_pops <- as.numeric(breseq_gene_props_day6$percent_pops)
breseq_gene_props_day6[breseq_gene_props_day6$percent_pops == 0, "percent_pops"] <- 10^-4

# Reemplazar guiones problemáticos y asegurar codificación adecuada
breseq_gene_props_day6$top_gene <- gsub("‑", "-", breseq_gene_props_day6$top_gene)  

# Mantener el orden de los factores en el eje X
breseq_gene_props_day6$top_gene <- factor(breseq_gene_props_day6$top_gene, levels = rev(unique(breseq_gene_props_day6$top_gene)))

#--- Graficar Fig 4B---- 
ggplot(breseq_gene_props_day6, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")  # Definir todos los ajustes en una sola llamada
  ) +
  scale_fill_viridis(
    name = "Poblaciones\ncon mutación (%)",
    limits = c(0, 100),  # Definir límites de escala
    breaks = c(0, 25, 50, 75, 100)  # Definir breaks específicos a mostrar
  ) +
  ylab("") +
  xlab("") +
  ggtitle("Día 6 de evolución experimental")

##--- Randomization Test (Day 36)----
# Construcción de la estadística de la prueba
tmp2_day36 <- sim_coef_NS_noP0[sim_coef_NS_noP0$passage1 == "P12" & sim_coef_NS_noP0$passage2 == "P12" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2), ]

test_same_P12_day36 <- tmp2_day36[tmp2_day36$same_gene == "same", "sim_coef"]
test_diff_P12_day36 <- tmp2_day36[tmp2_day36$same_gene == "different", "sim_coef"]
test_gene_P12_day36 <- mean(test_same_P12_day36) - mean(test_diff_P12_day36)

# Aleatorización y remuestreo
perm_gene_P12_day36 <- c()
set.seed(123)
for (i in seq(1, 10000)) {
  tmp2_day36$same_gene <- sample(tmp2_day36$same_gene, replace = FALSE)
  perm_same_P12_day36 <- mean(tmp2_day36[tmp2_day36$same_gene == "same", "sim_coef"])
  perm_diff_P12_day36 <- mean(tmp2_day36[tmp2_day36$same_gene == "different", "sim_coef"])
  perm_gene_P12_day36 <- c(perm_gene_P12_day36, perm_same_P12_day36 - perm_diff_P12_day36)
}

# P-valor (proporción de permutaciones más extremas que el valor observado)
p_value_P12_day36 <- length(perm_gene_P12_day36[perm_gene_P12_day36 >= test_gene_P12_day36]) / length(perm_gene_P12_day36)
print(p_value_P12_day36)

# Preparando datos para visualización
test_gene_df_day36 <- data.frame(
  sim_coef = c(test_same_P12_day36, test_diff_P12_day36),
  label = c(rep("Same\nresistance\ngene", length(test_same_P12_day36)), rep("Different\nresistance\ngene", length(test_diff_P12_day36)))
)
colnames(test_gene_df_day36) <- c("sim_coef", "label")
test_gene_df_day36$label <- factor(test_gene_df_day36$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene"))

#--- Fig 5A ---- 
# Creando y visualizando la gráfica sin anotación del p-valor
p_day36 <- ggplot(test_gene_df_day36, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Overlap in acquired mutations\n(Paired Sørensen-Dice similarity)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55)

print(p_day36)

## Identificando qué genes mutaron en qué poblaciones
# Asumiendo que 'breseq_NS_noP0' está cargado en el entorno
top_genes_in_study2_day36 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(top_genes_in_study2_day36) <- c("description", "Line", "hits")
top_genes_day36 <- names(tail(sort(table(top_genes_in_study2_day36[top_genes_in_study2_day36$hits > 0, "description"])), 20))

### Construyendo una tabla de análisis de genes para P12
breseq_gene_analysis2_day36 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(breseq_gene_analysis2_day36) <- c("gene", "Line", "hits")

# Definiendo los vectores de genes de resistencia
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  

# Anotar cada población según su gen de resistencia
breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% glycoside, "res_gene"] <- "glycoside"

# Eliminar filas donde res_gene es NA
breseq_gene_analysis2_day36 <- breseq_gene_analysis2_day36[!is.na(breseq_gene_analysis2_day36$res_gene), ]

### Calcular el porcentaje de poblaciones con cada gen de resistencia que mutó en cada uno de los otros genes
breseq_gene_props_day36 <- data.frame(matrix(nrow = 0, ncol = 3))
l_day36 <- list(rfbA = rfbA, glycosyl = glycosyl, glycoside = glycoside)
names(l_day36) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants")

# Bucle a través de genes y cálculo de proporciones de mutación
for (gene in top_genes_day36) {
  for (j in seq_along(l_day36)) {
    res_gene <- l_day36[[j]]
    prop <- nrow(breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% res_gene & 
                                               breseq_gene_analysis2_day36$gene == gene & 
                                               breseq_gene_analysis2_day36$hits > 0, ]) / 
      nrow(breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% res_gene & 
                                         breseq_gene_analysis2_day36$gene == gene, ])
    breseq_gene_props_day36 <- rbind(breseq_gene_props_day36, c(gene, names(l_day36)[j], prop))
  }
}

# Ordenar la columna 'res_gene' para que "rfbA mutants" aparezca primero, seguido de "PSPTO_4988 mutants" y luego "PSPTO_4991 mutants"
colnames(breseq_gene_props_day36)
colnames(breseq_gene_props_day36) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day36$res_gene <- factor(breseq_gene_props_day36$res_gene, levels = c("PSPTO_4991 mutants", "PSPTO_4988 mutants", "rfbA mutants"))

# Asignar nombres de columnas
colnames(breseq_gene_props_day36) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day36$percent_pops <- as.numeric(breseq_gene_props_day36$percent_pops)
breseq_gene_props_day36[breseq_gene_props_day36$percent_pops == 0, "percent_pops"] <- 10^-4

# Reemplazar guiones problemáticos y asegurar una codificación adecuada
breseq_gene_props_day36$top_gene <- gsub("‑", "-", breseq_gene_props_day36$top_gene)

# Mantener el orden de los factores en el eje X
breseq_gene_props_day36$top_gene <- factor(breseq_gene_props_day36$top_gene, levels = rev(unique(breseq_gene_props_day36$top_gene)))

#--- Gráfica Fig 5B----
ggplot(breseq_gene_props_day36, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")) +
  scale_fill_viridis(name = "Poblaciones\ncon mutación (%)", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  ylab("") +
  xlab("") +
  ggtitle("Día 36 de evolución experimental")

# Guardando objetos necesarios para las figuras
save(p_day6, p_day36, annotation_day6, test_gene_df_day6, annotation_day6, breseq_gene_props_day6, test_gene_df_day36, breseq_gene_props_day36, 
     file = "ReproHack_data_figures4_5.RData")

sessionInfo()
#R version 4.4.0 (2024-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 22631)

#Matrix products: default


#locale:
#  [1] LC_COLLATE=Spanish_Mexico.utf8  LC_CTYPE=Spanish_Mexico.utf8    LC_MONETARY=Spanish_Mexico.utf8
#[4] LC_NUMERIC=C                    LC_TIME=Spanish_Mexico.utf8    

#time zone: America/Mexico_City
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] viridis_0.6.5     viridisLite_0.4.2 ggplot2_3.5.1     readxl_1.4.3     

#loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       generics_0.1.3    labeling_0.4.3    glue_1.7.0       
#[7] colorspace_2.1-0  gridExtra_2.3     scales_1.3.0      fansi_1.0.6       grid_4.4.0        cellranger_1.1.0 
#[13] munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4   compiler_4.4.0    dplyr_1.1.4       pkgconfig_2.0.3  
#[19] rstudioapi_0.16.0 farver_2.1.2      R6_2.5.1          tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0     
#[25] magrittr_2.0.3    tools_4.4.0       withr_3.0.1       gtable_0.3.5   
