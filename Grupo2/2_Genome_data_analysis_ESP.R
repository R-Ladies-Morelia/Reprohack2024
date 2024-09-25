# Title: Analysis of bacterial genomes evolved in the absence of phage (# Análisis de genomas bacterianos evolucionados en ausencia de fagos)
# Author: Reena Debray
# Date: Feb 1, 2022
# Modified: Amairani Cancino Bello (04/Sep/2024)


###--- Functions ----
## Esta función utiliza la salida del programa breseq para construir una matriz de genes que fueron mutados en cada población durante la evolución experimental
### Toma datos en la siguiente forma: "breseq", un data frame con una entrada para cada mutación en cada población en cada momento y su frecuencia alélica correspondiente
### "gene_list", una lista de todos los genes observados en la población o subpoblación de interés
### "sample_list", una lista de todas las muestras en la población o subpoblación de interés 
### AF_min, una frecuencia alélica mínima para consideración. Establezco AF_min=0 (sin mínimo), aunque note que breseq solo devuelve polimorfismos con una frecuencia poblacional de 0.05 o mayor.


# --- Packages ----
#' @import readxl : This package is used to read Excel files (.xls and .xlsx).
#' @import ggplot2: This is one of the most popular packages in R for data visualization. You can create basic plots (such as histograms and bar charts) as well as more complex visualizations (like heatmaps and violin plots).
#' @import viridis : This package provides color palettes for plots in R. The viridis palettes are particularly useful in visualizations like heatmaps or scatter plots where color differentiation is important.

library(readxl)
BiocManager::install("ggplot2")
library(ggplot2)
install.packages("viridis")
library(viridis)


#--- Specify the full directory path ----
#' @title Especificar la Ruta del Directorio
#' @description Establece la ruta del directorio que contiene las tablas de mutación.
#' @param breseq Una cadena que especifica la ruta del directorio.
#' @usage breseq <- "C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables"

#Directorio que contiene las tablas de mutación.
breseq <- "C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables"

#--- Generate Annotared Gene Matrix---- 
#' @title Generar Matriz de Genes Anotada
#' @description Crea una matriz anotada con genes y muestras utilizando los datos de mutación breseq.
#' @param breseq Data frame que contiene los datos de mutación breseq.
#' @param gene_list Lista de genes a incluir en la matriz.
#' @param sample_list Lista de muestras a incluir en la matriz.
#' @param AF_min Frecuencia alélica mínima para el filtrado.
#' @return Una matriz anotada con el número de sitios de mutación por gen por muestra.
#' @export
gene_matrix_annotated<-function(breseq,gene_list,sample_list,AF_min){
  LOR_gene_matrix<-matrix(ncol=length(gene_list),nrow=length(sample_list)) #Se crea una matriz llamada LOR_gene_matrix con el número de filas igual a la cantidad de muestras (sample_list) y el número de columnas igual a la cantidad de genes (gene_list).
  colnames(LOR_gene_matrix)=gene_list #Se asignan nombres de columna (genes) 
  rownames(LOR_gene_matrix)=sample_list #Se asignan nombres de fila (muestras).
  
  for (gene in gene_list){
    for (sample in sample_list){
      if (gene%in%breseq[breseq$Sample==sample & breseq$freq>=AF_min,"gene_2"]){num_sites<-length(unique(breseq[breseq$Sample==sample & breseq$gene_2==gene & !is.na(breseq$gene_2),"position"])); LOR_gene_matrix[sample,gene]=num_sites}
      else{LOR_gene_matrix[sample,gene]=0}
    }
  }
  return(LOR_gene_matrix)
} #La función recorre cada gen en la lista `gene_list` y, para cada uno, evalúa las mutaciones en todas las muestras de `sample_list`. 
#Verifica si el gen está presente en el `data frame` de mutaciones (`breseq`) y si su frecuencia alélica es mayor o igual al umbral `AF_min`. 
#Si el gen aparece en la muestra, cuenta las posiciones únicas mutadas y asigna este valor en la matriz correspondiente; de lo contrario, asigna un 0.


#---Read and filter breseq output----
#Leer multiples archivos de excel y almacenarlos en un solo DataFrame
#' @title Leer y Filtrar Salida de breseq
#' @description Lee los archivos de datos de mutación breseq, formatea el data frame y filtra los datos de mutación basándose en varios criterios.
#' @return Un data frame que contiene datos de mutación filtrados.
#' @export

breseq_annotated<-data.frame() # data frame vacío para almacenar la información combinada de los archivos.
filenames=list.files("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables") #se genera una lista con los nombres de todos los archivos en la carpeta especificada.
for (file in filenames){
  ann <- data.frame(read_excel(paste("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables", file, sep="/")))
  ann$Line<-unlist(strsplit(file,split="_"))[1]  #Itera sobre cada archivo, lee cada archivo de Excel y lo convierte en un data frame.
  ann$Passage<-unlist(strsplit(file,split="_"))[2]
  ann$Sample<-paste(ann$Line,ann$Passage,sep="_")
  if (ncol(ann)>10){ann<-ann[,-c(8:12)]} # some files have an extra column of notes for annotations without evidence (*) resulting in blanks for cols 8-12
  breseq_annotated<-rbind(breseq_annotated,ann) #Se usa rbind para añadir los datos del archivo actual a breseq_annotated.
} 
dim(breseq_annotated) #DataFrame resultante 

#--- Remove unassigned evidence----
#' @title Eliminar Evidencia No Asignada
#' @description Filtra la evidencia no asignada de los datos de mutación breseq.
#' @param breseq_annotated Data frame que contiene datos de mutación.
#' @return Data frame filtrado sin evidencia no asignada.
#' @export

breseq_annotated<-breseq_annotated[!breseq_annotated$evidence%in%c(NA,"Unassigned missing coverage evidence","*","Unassigned new junction evidence","?"),]
dim(breseq_annotated)

# --- Function: Remove Sites Differing from Reference ----
#' @title Eliminar Sitios que Difieren del Referencia
#' @description Elimina cualquier sitio de mutación que difiera de la referencia en la línea ANCDC3000 (ancestral DC3000).
#' @param breseq_annotated Data frame. Contiene los datos anotados de breseq.
#' @return Data frame sin sitios que difieran de la referencia en ANCDC3000 (ancestral DC3000).

## Eliminar cualquier sitio que difiera de la referencia en ANCDC3000 (ancestral DC3000)
ANCDC3000_sites<-breseq_annotated[breseq_annotated$Line=="ANCDC3000","position"] #Se extraen todas las posiciones (position) donde la línea es "ANCDC3000"
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%ANCDC3000_sites,] #Se filtra el data frame breseq_annotated para eliminar cualquier fila donde la columna position coincida con alguna de las posiciones almacenadas en ANCDC3000_sites.
dim(breseq_annotated)


# --- Function: Remove Fixed Sites ----
#' @title Eliminar Sitios Fijos
#' @description Elimina sitios que están fijos en cada línea.
#' @param breseq_annotated Data frame. Contiene los datos anotados de breseq.
#' @return Data frame sin sitios de mutación fijos.

## Eliminar sitios que están fijos en cada línea
tab<-data.frame(table(breseq_annotated[breseq_annotated$freq==1,"Line"],breseq_annotated[breseq_annotated$freq==1,"position"])) # tabla de sitios fijos, Freq =1. Esto significa que la mutación está fija en esa línea.
tab<-tab[tab$Freq>0,] #Filtra la tabla tab para eliminar las combinaciones línea-mutación que no tienen ninguna mutación (Freq > 0).
mut_tab<-data.frame(table(tab$Var2))  #Se crea otra tabla que cuenta en cuántas líneas aparece cada mutación. La columna Var2 representa las posiciones de mutación.
sites_to_remove<-as.character(mut_tab[mut_tab$Freq==length(unique(breseq_annotated$Line)),"Var1"])
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%sites_to_remove,] #Filtra el data frame original breseq_annotated para eliminar las posiciones de mutación que están en el vector -sites_to_remove-, es decir, los sitios fijos.
dim(breseq_annotated)

# --- Function: Restructure Frameshift Mutations ----
#' @title Restructurar Mutaciones por Desplazamiento de Lectura
#' @description Restructura las mutaciones por desplazamiento de lectura para ser consistentes con las mutaciones puntuales.
#' @param breseq_annotated Data frame. Contiene los datos anotados de breseq.
#' @return Data frame con las mutaciones por desplazamiento de lectura restructuradas.

## Restructurar las mutaciones por desplazamiento de lectura para ser consistentes con las mutaciones puntuales
for (i in seq(1,nrow(breseq_annotated))){ 
  pos<-breseq_annotated[i,"position"] 
  if (grepl(":",pos,fixed=TRUE)){
    pos<-unlist(strsplit(pos,split=":"))[1]
    pos<-paste(unlist(strsplit(pos,split=",")),collapse="")
    breseq_annotated[i,"position"]<-pos
  }
}
breseq_annotated$position<-as.numeric(breseq_annotated$position) 


# --- Function: Restructure Gene Column ----
#' @title Restructurar Columna de Genes
#' @description Restructura la columna de genes para eliminar la dirección.
#' @param breseq_annotated Data frame. Contiene los datos anotados de breseq.
#' @return Data frame con la columna de genes restructurada.

## Restructurar la columna de genes para eliminar la dirección
for (i in seq(1,nrow(breseq_annotated))){
  breseq_annotated[i,"gene_2"]<-substr(breseq_annotated[i,"gene"],1,nchar(breseq_annotated[i,"gene"])-2) # recorre cada fila y elimina los últimos dos caracteres
}

# --- Function: Filter Mutations with Sliding Window ----
#' @title Filtrar Mutaciones con Ventana Deslizante
#' @description Filtra las mutaciones utilizando una ventana deslizante de 50 bp dentro de cada línea/punto de tiempo.
#' @param breseq_annotated Data frame. Contiene los datos anotados de breseq.
#' @param N Numérico. Número máximo de vecinos permitidos dentro de la ventana deslizante.
#' @return Data frame con las mutaciones filtradas.

#Si una ventana contiene más de N mutaciones, todas las mutaciones en esa ventana se eliminan.
breseq_annotated_filtered<-data.frame() #data frame vacío 
N=4  #establece que si hay más de 4 mutaciones dentro de una ventana de 50 bp, las mutaciones serán eliminadas.


for (sample in unique(breseq_annotated$Sample)){
  pre_filter_POS<-sort(as.numeric(breseq_annotated[breseq_annotated$Sample==sample,"position"])) #se ordenan las posiciones de las mutaciones (pre_filter_POS) en formato numérico.
  sites_to_remove<-c() # -sites_to_remove- almacenará las posiciones que deben ser eliminadas.
  for (pos in pre_filter_POS){
    # ventana deslizante
    min=pos-50
    while (min<pos){
      max<-min+50
      if (length(pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max])>N)   {sites_to_remove<-unique(c(sites_to_remove,pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max]))} # si una ventana contiene más de 3 mutaciones, elimine todas
      min<-min+1
    }
  }
  
  #' @title Filtrar y Anotar Datos Genéticos
  #'
  #' @Description Esta función filtra los sitios basándose en una muestra proporcionada y una lista de sitios a eliminar, y luego categoriza las anotaciones en los datos filtrados.
  #'
  #' @param breseq_annotated Un data frame que contiene los datos genéticos anotados.
  #' @param sample Una cadena que indica la muestra por la que filtrar.
  #' @param sites_to_remove Un vector de posiciones a eliminar de los datos.
  #'
  #' @return Un data frame con anotaciones filtradas y tipos actualizados.
  
  ### filtrar sitios
  post_filter<-breseq_annotated[breseq_annotated$Sample==sample & !breseq_annotated$position%in%sites_to_remove,] #seleccionar únicamente las filas correspondientes a la muestra indicada por sample y eliminar aquellas posiciones de sitios que estén en la lista sites_to_remove.
  breseq_annotated_filtered<-rbind(breseq_annotated_filtered,post_filter) # data frame resultante (post_filter) se añade al data frame global breseq_annotated_filtered utilizando la función rbind(), lo que permite acumular los datos filtrados para cada muestra en un solo objeto.
}
# Imprimir las dimensiones de los datos filtrados
dim(breseq_annotated_filtered)

## Resumir anotaciones como NS, S, pseudogénico o intergénico
for (i in seq(1,nrow(breseq_annotated_filtered))){
  anno<-substr(breseq_annotated_filtered[i,"annotation"],1,5) #Extraer los primeros 5 caracteres de la columna "annotation"
  if (anno=="pseud"){breseq_annotated_filtered[i,"Type"]="Pseudogene"} # Si la anotación comienza con "pseud", se clasifica como pseudogénico
  else if (anno=="inter"){breseq_annotated_filtered[i,"Type"]="Intergenic"} #Si la anotación comienza con "inter", se clasifica como intergénico
  else { # Extraer el primer y último carácter de la cadena de 5 caracteres (representan aminoácidos)
    AA1<-substr(anno,1,1)  # Primer aminoácido 
    AA2<-substr(anno,5,5)  # Último aminoácido
    if (AA1==AA2){breseq_annotated_filtered[i,"Type"]="S"} # Si el primer y último aminoácido son iguales, la mutación es sinónima (S)
    else {breseq_annotated_filtered[i,"Type"]="NS"}   # Si son diferentes, la mutación es no sinónima (NS)
  }
}  

## Contar número de sitios fijos y polimorfismos en cada población en el pase 12

#' @title Contar Sitios Fijos y Polimorfismos en Cada Población en el Pase 12
#'
#' @description Este script calcula el número de mutaciones fijas y polimorfismos en cada población en el pase 12, excluyendo las mutaciones de resistencia adquiridas antes de la evolución experimental. También realiza análisis estadísticos para comparar los conteos entre diferentes poblaciones e identifica poblaciones mutadoras.
#'
#' @details
#' El script incluye los siguientes análisis:
#' \itemize{
#'   \item Cuenta el número de mutaciones fijas (tanto sinónimas como no sinónimas) en cada población en el pase 12, excluyendo las mutaciones de resistencia.
#'   \item Calcula la media y la desviación estándar del conteo de mutaciones fijas y realiza una prueba t para comparar el número de mutaciones fijas entre diferentes poblaciones.
#'   \item Cuenta el número de polimorfismos (tanto sinónimos como no sinónimos) en cada población en el pase 12, excluyendo las mutaciones de resistencia.
#'   \item Calcula la media y la desviación estándar del conteo de polimorfismos y realiza una prueba t para comparar el número de polimorfismos entre diferentes poblaciones.
#'   \item Identifica poblaciones mutadoras basándose en criterios específicos y realiza una prueba t para comparar el número de mutaciones fijas entre poblaciones mutadoras y no mutadoras.
#' }
#'
#' @importFrom readxl read_excel
#' @param Costs_of_res Un data frame que contiene posiciones de mutaciones de resistencia a excluir del conteo.
#' @param breseq_annotated_filtered Un data frame que contiene datos de mutaciones, incluyendo información de pase, tipo, posición, frecuencia y línea.
#' @return Imprime la media, la desviación estándar y los resultados de las pruebas t para el número de mutaciones fijas y polimorfismos en cada población, y proporciona comparaciones entre poblaciones mutadoras y no mutadoras.
#' @examples
#' \dontrun{
#' # Cargar los datos necesarios y ejecutar el script para calcular los conteos y realizar los análisis
#' }
#' @export

### Mutaciones fijas (todos los tipos)
### Excluir mutaciones de resistencia del conteo (porque fueron adquiridas antes de la evolución experimental)
Costs_of_res <- read_excel("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Costs_of_res.xlsx") #Cargar el archivo de Excel con las mutaciones de resistencia.
res_mutations<-unique(Costs_of_res$Position) # Extrae las posiciones de las mutaciones de resistencia y las almacena en un vector -res_mutations-.
#Filtrar las mutaciones que ocurrieron en el pase 12, filtrar las mutaciones que son sinónimas (S) o no sinónimas (NS) y excluir las mutaciones de resistencia (aquellas que están en res_mutations).
a<-aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"freq"],by=list(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"Line"]),FUN=function(x){length(x[x==1])})
mean(a$x) #Calcula la media del número de mutaciones fijas por población.
sd(a$x) #Calcula la desviación estándar del número de mutaciones fijas por población.
t.test(x~(substr(a$Group.1,1,3)=="ANC"),a,var.equal=T) # Realiza una prueba t para comparar el número de mutaciones fijas entre diferentes poblaciones.

### Polimorfismos (todos los tipos)
p<-aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"freq"],by=list(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"Line"]),FUN=function(x){length(x[x<1])})
mean(p$x)
sd(p$x) 
t.test(x~(substr(p$Group.1,1,3)=="ANC"),p,var.equal=T)

### Mutaciones fijas en poblaciones mutadoras (dentro de genes)
mut<-c("FMS13","MS1","MS10","MS15","QAC4") #Lista que contiene los nombres de las poblaciones que han sido identificadas como mutantes.
p[p$Group.1%in%mut,"Mutador"]<-"Y" #Asigna "Y" a las poblaciones que están en la lista mut y "N" a las que no lo están.
p[!p$Group.1%in%mut,"Mutador"]<-"N"
t.test(x~Mutador,p,alternative="less",var.equal=T) #Realizar una prueba t para comparar el número de polimorfismos entre las poblaciones (Y) y (N).

#--- Figure 4: Genomic parallelism across populations---- 

#' @title Paralelismo Genómico Entre Poblaciones: Análisis de Similaridad Sørensen-Dice
#'
#' @description Este script analiza el paralelismo genómico entre poblaciones iniciadas con bacterias resistentes a fagos. Identifica mutaciones no sinónimas que no están fijas en el pase 0, construye una matriz de genes, calcula los coeficientes de similitud Sørensen-Dice para todos los pares de muestras y anota los genes de resistencia de las poblaciones.
#'
#' @details
#' El script realiza los siguientes pasos:
#' \itemize{
#'   \item Identifica poblaciones iniciadas con bacterias resistentes a fagos, excluyendo controles ancestrales.
#'   \item Construye un data frame de mutaciones no sinónimas que no están fijas en el pase 0 en estas poblaciones.
#'   \item Crea una matriz de genes utilizando solo mutaciones no fijas en el pase 0.
#'   \item Calcula los coeficientes de similitud Sørensen-Dice para cada par de muestras para evaluar el grado de superposición en las mutaciones genéticas entre muestras.
#'   \item Anota los genes de resistencia de las poblaciones y determina si los pares de muestras tienen los mismos o diferentes genes de resistencia.
#' }
#'
#' @param breseq_annotated_filtered Un data frame que contiene datos anotados de mutaciones, incluyendo línea, pase, posición, tipo, frecuencia, gen y información de muestra.
#' @param Costs_of_res Un data frame que contiene información sobre los genes de resistencia asociados con diferentes poblaciones.
#' @return Un data frame (`sim_coef_NS_noP0`) que contiene los coeficientes de similitud Sørensen-Dice para cada par de muestras, anotado con información de genes de resistencia y si las muestras comparten los mismos o diferentes genes de resistencia.
#' @examples
#' \dontrun{
#' # Ejecutar el script para calcular los coeficientes de similitud y anotar los genes de resistencia entre poblaciones
#' }
#' @export


### Generate names of populations initiated with phage-resistant bacteria only (exclude ancestral controls)
n<-unique(breseq_annotated_filtered[substr(breseq_annotated_filtered$Line,1,3)!="ANC","Line"]) # Generando un vector "n" con los nombres de las poblaciones de interés que no son ancestrales.

### Form a data frame of non-synonymous mutations not fixed at Psg 0 only
breseq_NS_noP0<-data.frame(matrix(nrow=0,ncol=ncol(breseq_annotated_filtered))) #Crea un data frame vacío con el mismo número de columnas que el data frame original
for (line in n){ #Itera sobre cada población identificada.
  P0_mutations<-breseq_annotated_filtered[breseq_annotated_filtered$Line==line & breseq_annotated_filtered$Passage=="P0" & breseq_annotated_filtered$freq==1,"position"] #Obtiene las posiciones de mutaciones que están fijas (frecuencia igual a 1) en el pase 0.
  ###rows corresponding to P2/P12 mutations not fixed in P0 of this sample
  sample_data<-breseq_annotated_filtered[breseq_annotated_filtered$Line==line & breseq_annotated_filtered$Passage!="P0" & !breseq_annotated_filtered$position%in%P0_mutations & breseq_annotated_filtered$Type=="NS",] #Selecciona las mutaciones no sinónimas (Type == "NS") que no están fijas en el pase 0
  
  breseq_NS_noP0<-rbind(breseq_NS_noP0,sample_data) # Agregar las mutaciones relevantes a breseq_NS_noP0.
}

### Make a gene matrix using only mutations not fixed in P0
gene_list<-sort(unique(breseq_NS_noP0$"gene_2")) #Extraer la columna gene_2 del data frame breseq_NS_noP0, que contiene los nombres de los genes asociados con las mutaciones. unique= Obtiene los nombres únicos de los genes, eliminando duplicados.sort=Ordena alfabéticamente la lista de genes.
sample_list<-unique(breseq_NS_noP0$Sample) #Extrae la columna Sample, que contiene las identificaciones de las muestras.
LOR_gene_matrix_NS_noP0<-gene_matrix_annotated(breseq_NS_noP0,gene_list,sample_list,0) 

### Calculate Sørensen-Dice similarity coefficients for every pair of samples
all_pairs<-combn(sample_list,2) #Genera todas las combinaciones posibles de pares de muestras (2 a la vez).
sim_coef_NS_noP0<-data.frame(matrix(nrow=0,ncol=3)) #Crea un data frame vacío con 3 columnas y sin filas. Este data frame se llenará con los resultados de similitud.

for (i in seq(1,ncol(all_pairs))){ #Itera sobre cada columna en all_pairs, que contiene los pares de muestras.
  sample1<-all_pairs[1,i] #Asigna las dos muestras del par actual.
  sample2<-all_pairs[2,i]
  
  ### Names of genes with any mutations from ancestral genotype (any number of independent mutations per gene)
  #obtienen los genes con mutaciones y se verifica qué genes tienen mutaciones (valores > 0).
  sample1_genes<-names(LOR_gene_matrix_NS_noP0[sample1,apply(data.frame(LOR_gene_matrix_NS_noP0[sample1,]),1,function(x){x>0})])
  sample2_genes<-names(LOR_gene_matrix_NS_noP0[sample2,apply(data.frame(LOR_gene_matrix_NS_noP0[sample2,]),1,function(x){x>0})])
  
  ### Calculate the extent of overlap
  #Calcular la extensión del solapamiento
  sim_coef<-(2*length(intersect(sample1_genes,sample2_genes))) / (length(sample1_genes) + length(sample2_genes))
  sim_coef_NS_noP0<-rbind(sim_coef_NS_noP0,c(sample1,sample2,sim_coef))
}
colnames(sim_coef_NS_noP0)=c("sample1","sample2","sim_coef")
sim_coef_NS_noP0$sim_coef<-as.numeric(sim_coef_NS_noP0$sim_coef) #sim_coef almacena el valor de similitud Sørensen-Dice, que varía entre 0 (sin superposición) y 1 (muestras idénticas en términos de genes mutados).

### Annotate the resistance genes of the populations
#Agregan cuatro nuevas columnas al data frame sim_coef_NS_noP0:
#line1 y passage1 para la primera muestra.
#line2 y passage2 para la segunda muestra.

for (i in seq(1,nrow(sim_coef_NS_noP0))){
  sim_coef_NS_noP0[i,"line1"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample1"],split="_"))[1]
  sim_coef_NS_noP0[i,"passage1"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample1"],split="_"))[2]
  sim_coef_NS_noP0[i,"line2"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample2"],split="_"))[1]
  sim_coef_NS_noP0[i,"passage2"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample2"],split="_"))[2]
 # buscar en el data frame Costs_of_res la fila correspondiente a la población específica (definida en line1) y extrae el nombre del gen de resistencia asociado. Se hace lo mismo con line2 para la segunda muestra.
  sim_coef_NS_noP0[i,"gene1"]=as.character(unique(Costs_of_res[Costs_of_res$Population==sim_coef_NS_noP0[i,"line1"],"Gene"]))
  sim_coef_NS_noP0[i,"gene2"]=as.character(unique(Costs_of_res[Costs_of_res$Population==sim_coef_NS_noP0[i,"line2"],"Gene"]))
}

# Annotate whether pairs have the same or different resistance genes
#Verificar que tanto gene1 como gene2 no sean valores NA (es decir, que haya un gen anotado para ambas muestras).Comparar si los genes de resistencia de las dos muestras son iguales.
#Si las dos condiciones se cumplen, asigna el valor "same" a la nueva columna same_gene. Si los genes no son iguales se agrega "different".
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1==sim_coef_NS_noP0$gene2,"same_gene"]="same"
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1!=sim_coef_NS_noP0$gene2,"same_gene"]="different"

##--- Randomization test (Day 6)---- 
#' @title Prueba de Aleatorización y Visualización para el Análisis de Paralelismo Genómico (Día 6)
#'
#' @description Este script realiza una prueba de aleatorización para evaluar la significancia de la superposición en mutaciones adquiridas entre poblaciones, basada en los coeficientes de similitud Sørensen-Dice. Construye una estadística de prueba, realiza aleatorización y re-muestreo, calcula el valor p y visualiza los resultados utilizando un diagrama de caja.
#'
#' @details 
#' El script realiza los siguientes pasos:
#' \itemize{
#'   \item Construye una estadística de prueba que representa la diferencia en los coeficientes de similitud promedio entre pares de muestras con el mismo gen de resistencia y aquellos con genes de resistencia diferentes.
#'   \item Realiza una prueba de aleatorización permutando las etiquetas ("mismo" o "diferente" gen de resistencia) y recalcula la estadística de prueba 10,000 veces para generar una distribución de la estadística de prueba bajo la hipótesis nula.
#'   \item Calcula el valor p como la proporción de estadísticas de prueba permutadas que son más extremas que la estadística de prueba observada.
#'   \item Prepara los datos para la visualización y crea un diagrama de caja para ilustrar la superposición en mutaciones adquiridas entre muestras con el mismo o diferente gen de resistencia.
#'   \item Anota el gráfico con el valor p calculado.
#' }
#'
#' @param sim_coef_NS_noP0 Un data frame que contiene los coeficientes de similitud Sørensen-Dice para cada par de muestras, anotado con información sobre el gen de resistencia.
#' @return Un gráfico que muestra los coeficientes de similitud Sørensen-Dice para pares de muestras con el mismo o diferente gen de resistencia, y el valor p para la prueba de aleatorización.
#' @examples
#' \dontrun{
#' # Ejecutar el script para realizar la prueba de aleatorización y generar el gráfico
#' }
#' @export


### Construct test statistic
# Filtrar el data frame sim_coef_NS_noP0 para analizar únicamente las muestras del pase 2 (P2), excluyendo aquellas donde los genes de resistencia (gene1 y gene2) sean NA.
tmp<-sim_coef_NS_noP0[sim_coef_NS_noP0$passage1=="P2" & sim_coef_NS_noP0$passage2=="P2" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2),]
test_same_P2<-tmp[tmp$same_gene=="same","sim_coef"] #Extraer los coeficientes de similitud Sørensen-Dice para los pares de muestras que comparten el mismo gen de resistencia (same_gene == "same").
test_diff_P2<-tmp[tmp$same_gene=="different","sim_coef"] #Extraer los coeficientes de similitud para los pares de muestras que tienen genes de resistencia diferentes (same_gene == "different").
test_gene_P2<-mean(test_same_P2)-mean(test_diff_P2) #promedio en los coeficientes de similitud entre los pares de muestras con el mismo gen de resistencia y aquellos con genes de resistencia diferentes.

### Randomization & resampling
perm_gene_P2<-c() # Crear un vector vacío para almacenar los resultados de las permutaciones.
set.seed(123) # Fijar la semilla para la reproducibilidad.
for (i in seq(1,10000)){ # Realizar 10,000 permutaciones.
  tmp$same_gene<-sample(tmp$same_gene,replace=F)  # Barajar la columna same_gene sin reemplazo 
  perm_same_P2<-mean(tmp[tmp$same_gene=="same","sim_coef"]) # Calcular el promedio de los coeficientes de similitud para los pares con el mismo gen de resistencia.
  perm_diff_P2<-mean(tmp[tmp$same_gene=="different","sim_coef"]) # Calcular el promedio de los coeficientes de similitud para los pares con genes de resistencia diferentes.
  perm_gene_P2<-c(perm_gene_P2,perm_same_P2-perm_diff_P2) # Almacenar la diferencia en los coeficientes de similitud de esta permutación.
}

### p-value (proportion of permutations more extreme than observed value)
print(length(perm_gene_P2[perm_gene_P2>=test_gene_P2])/length(perm_gene_P2)) # Calcular y mostrar el valor p.

#--- Modified to add value p----

# Prepare data
# Crear un data frame con los coeficientes de similitud y etiquetas para las categorías "mismo gen de resistencia" y "diferente gen de resistencia".
test_gene_df <- data.frame(
  sim_coef = c(test_same_P2, test_diff_P2),
  label = c(rep("Same\nresistance\ngene", length(test_same_P2)), rep("Different\nresistance\ngene", length(test_diff_P2)))
)
colnames(test_gene_df) <- c("sim_coef", "label") # Asignar nombres de columna apropiados. 
test_gene_df$label <- factor(test_gene_df$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene")) # Convertir las etiquetas a factor, especificando el orden de los niveles.

# Calculate day 6 p-value
p_value_P2 <- length(perm_gene_P2[perm_gene_P2 >= test_gene_P2]) / length(perm_gene_P2) # Calcular el valor p para el día 6.
print(p_value_P2) # Mostrar el valor p.

# Set up p-value annotation
#Crear una etiqueta de anotación que incluya el valor p en el gráfico. Si el valor p es extremadamente bajo
#(por debajo de 0.001), se anota como "p < 0.001", de lo contrario se muestra el valor exacto.
if (p_value_P2 < 0.001) {
  annotation <- "p < 0.001***"
} else {
  annotation <- paste0("p = ", format(p_value_P2, digits = 3))
}

#--- Fig. 4A ---- 
# Creation and display of the graph
p <- ggplot(test_gene_df, aes(label, sim_coef)) + # Creación de un gráfico utilizando ggplot, donde "label" se mapea en el eje x y "sim_coef" en el eje y.
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) + # Añadir un boxplot sin mostrar valores atípicos, con un ancho de caja ajustado y un grosor de línea de 0.8.
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) + # Añadir puntos dispersos (jitter) con un desplazamiento horizontal pequeño (width = 0.05), tamaño 2 y transparencia 0.2.
  theme_classic(base_size = 18) +  # Aplicar el tema clásico para el gráfico con un tamaño de texto base de 18.
  xlab("") +   # Elimina la etiqueta del eje x.
  ylab("Overlap in acquired mutations\n(Pairwise Sørensen-Dice similarity)") + # Añadir la etiqueta al eje y, describiendo la similitud en mutaciones adquiridas (coeficiente de Sørensen-Dice).

  theme(panel.spacing = unit(3, "lines"), 
        strip.text = element_text(size = 18, face = "bold")) + # Ajusta el espacio entre paneles y el tamaño y estilo del texto de las etiquetas del gráfico.
  ylim(0.1, 0.55) + # Establece los límites del eje y entre 0.1 y 0.55.
  annotate("text", x = 1.5, y = 0.55, label = annotation, size = 6, hjust = 0.5, vjust = 1.5, color = "black")  # Añade una anotación de texto en la posición x = 1.5 y = 0.55 con el valor p, ajusta el tamaño y la alineación del texto.

print(p) # Muestra el gráfico.
 
#--- Fig.4B ---- 
#' @title Fig. 4B: Mapa de Calor de Porcentajes de Mutaciones por Gen de Resistencia
#'
#' @description Este script genera un mapa de calor que muestra el porcentaje de poblaciones con cada gen de resistencia que tienen mutaciones en genes específicos principales para el Día 6 de evolución experimental. 
#' El mapa de calor visualiza las proporciones de poblaciones con mutaciones en los 20 genes principales, categorizados por su gen de resistencia.
#'
#' @details 
#' El script realiza los siguientes pasos:
#' \itemize{
#'   \item Identifica los 20 genes principales mutados en las poblaciones al Día 6 de evolución experimental.
#'   \item Construye un data frame de mutaciones que aparecieron para el Día 6, excluyendo mutaciones fijas antes del inicio del experimento y mutaciones de resistencia a fagos.
#'   \item Anota cada población según su gen de resistencia: "rfbA", "glycosyl", o "glycoside".
#'   \item Calcula el porcentaje de poblaciones con cada gen de resistencia que tienen mutaciones en cada uno de los genes principales.
#'   \item Prepara y traza un mapa de calor mostrando el porcentaje de poblaciones con mutaciones en los genes principales, utilizando una escala de colores para representar el porcentaje de poblaciones.
#' }
#'
#' @param breseq_NS_noP0 Un data frame que contiene información sobre mutaciones no sinónimas no fijas en el pase 0, incluyendo descripciones de genes, líneas y conteos de mutaciones.
#' @param Costs_of_res Un data frame que contiene información sobre los genes de resistencia asociados con poblaciones específicas.
#' @return Un mapa de calor que muestra el porcentaje de poblaciones con cada gen de resistencia que tienen mutaciones en los genes principales para el Día 6.
#' @examples
#' \dontrun{
#' # Ejecutar el script para generar el mapa de calor para los 20 genes principales
#' }
#' @export

## Identificar qué genes fueron mutados en qué poblaciones
# Crear un data frame que cuente cuántas veces aparece cada descripción de gen por cada línea (Line) en el pase 2 (P2).
top_genes_in_study <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P2", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P2", "Line"]))
colnames(top_genes_in_study) = c("description", "Line", "hits") # Renombrar las columnas como "description", "Line", y "hits" (frecuencia de aparición).
top_genes <- names(tail(sort(table(top_genes_in_study[top_genes_in_study$hits > 0, "description"])), 20)) # Identificar los 20 genes más comunes en el estudio, excluyendo aquellos que no tienen hits.

### Construir un data frame de mutaciones que aparecieron para el Día 6, excluyendo mutaciones de resistencia a fagos o cualquier otra mutación fija antes del inicio del experimento
breseq_gene_analysis <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P2", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P2", "Line"]))
colnames(breseq_gene_analysis) = c("gene", "Line", "hits") # Renombrar las columnas como "gene", "Line", y "hits" para el nuevo data frame.

# Definir los vectores o listas de valores correspondientes a cada gen de resistencia
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")# Definir las líneas que tienen el gen de resistencia "rfbA"
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1") # Definir las líneas que tienen el gen de resistencia "glycosyl".
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  # Definir las líneas que tienen el gen de resistencia "glycoside".


### Anotar cada población por su gen de resistencia
breseq_gene_analysis[breseq_gene_analysis$Line %in% rfbA, "res_gene"] <- "rfbA" # Asignar "rfbA" como el gen de resistencia para las líneas que pertenecen a rfbA.
breseq_gene_analysis[breseq_gene_analysis$Line %in% glycosyl, "res_gene"] <- "glycosyl" # Asignar "glycosyl" como el gen de resistencia para las líneas que pertenecen a glycosyl.
breseq_gene_analysis[breseq_gene_analysis$Line %in% glycoside, "res_gene"] <- "glycoside" # Asignar "glycoside" como el gen de resistencia para las líneas que pertenecen a glycoside.

# Eliminar filas donde res_gene es NA
breseq_gene_analysis <- breseq_gene_analysis[!is.na(breseq_gene_analysis$res_gene), ] # Eliminar filas donde no se haya asignado un gen de resistencia (NA en res_gene).

### Calcular el porcentaje de poblaciones con cada gen de resistencia que mutó en cada uno de los otros genes
breseq_gene_props <- data.frame(matrix(nrow = 0, ncol = 3)) # Crear un data frame vacío con 3 columnas para almacenar los resultados: gen, gen de resistencia, y porcentaje de poblaciones.
l <- list(rfbA, glycosyl, glycoside) # Crear una lista con los grupos de líneas asociadas a cada gen de resistencia.
names(l) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants") # Asignar nombres a cada grupo de líneas de resistencia.

# Iterar a través de los genes y calcular las proporciones de mutaciones
for (gene in top_genes) { # Iterar sobre cada gen en la lista de los genes más comunes.
  for (j in seq_along(l)) { # Iterar sobre cada grupo de genes de resistencia (rfbA, glycosyl, glycoside).
    res_gene <- l[[j]] # Obtener las líneas correspondientes al gen de resistencia actual.
    #Calcular la proporción de poblaciones que tienen una mutación en el gen y en la línea correspondiente. 
    prop <- nrow(breseq_gene_analysis[breseq_gene_analysis$Line %in% res_gene &  
                                        breseq_gene_analysis$gene == gene & 
                                        breseq_gene_analysis$hits > 0, ]) / 
      nrow(breseq_gene_analysis[breseq_gene_analysis$Line %in% res_gene & 
                                  breseq_gene_analysis$gene == gene, ])
    breseq_gene_props <- rbind(breseq_gene_props, c(gene, names(l)[j], prop)) # Añadir el resultado al data frame breseq_gene_props.
  }
}

# Ordenar la columna 'res_gene' para que "rfbA mutants" aparezca primero, seguido de "PSPTO_4988 mutants" y luego "PSPTO_4991 mutants"
colnames(breseq_gene_props) # Mostrar los nombres actuales de las columnas.
colnames(breseq_gene_props) <- c("top_gene", "res_gene", "percent_pops") # Asignar los nombres correctos a las columnas del data frame.
# Convertir la columna 'res_gene' en un factor con un orden específico para que las categorías aparezcan en el orden deseado: "PSPTO_4991 mutants", "PSPTO_4988 mutants", y "rfbA mutants".
breseq_gene_props$res_gene <- factor(breseq_gene_props$res_gene, levels = c("PSPTO_4991 mutants", "PSPTO_4988 mutants", "rfbA mutants"))
breseq_gene_props$percent_pops <- as.numeric(breseq_gene_props$percent_pops) # Convertir la columna 'percent_pops' a tipo numérico para asegurar un formato adecuado para cálculos.
# Reemplazar cualquier valor de 0 en la columna 'percent_pops' con un valor muy pequeño (10^-4) para evitar problemas en las visualizaciones o cálculos futuros.
breseq_gene_props[breseq_gene_props$percent_pops == 0, "percent_pops"] <- 10^-4 

# Reemplazar guiones problemáticos y asegurar una codificación adecuada
breseq_gene_props$top_gene <- gsub("‑", "-", breseq_gene_props$top_gene)  

# Mantener el orden de los factores en el eje X
breseq_gene_props$top_gene <- factor(breseq_gene_props$top_gene, levels = rev(unique(breseq_gene_props$top_gene)))

#--- Plot Fig 4B---- 
# Crear un gráfico utilizando ggplot2, donde las variables 'top_gene' y 'res_gene' se usan en los ejes, 
#y el color de las celdas (tiles) se basa en 'percent_pops' multiplicado por 100 para expresarlo en porcentaje.
ggplot(breseq_gene_props, aes(top_gene, res_gene, fill = percent_pops * 100)) + 
  geom_tile(color = "grey80") + # Agregar celdas rectangulares (tiles) al gráfico, con los bordes en color gris claro.
  theme_classic(base_size = 14) + # Aplicar el tema clásico de ggplot, que elimina el fondo gris y simplifica el gráfico, con un tamaño base de fuente de 14.
  # Ajustar el texto del eje X para que esté en ángulo de 90 grados, con alineación horizontal a 1, tamaño de fuente 10 y la familia tipográfica "Arial".
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")  
  ) +
  scale_fill_viridis(
    name = "Populations\nwith mutation (%)", # Definir el nombre de la escala de color.
    limits = c(0, 100),# Establecer los límites de la escala entre 0 y 100%.
    breaks = c(0, 25, 50, 75, 100)  # Definir los puntos específicos en los que aparecerán las etiquetas de la escala.
  ) +
  ylab("") + 
  xlab("") +
  ggtitle("Day 6 of experimental evolution")  # Agregar el título "Day 6 of experimental evolution" al gráfico.


##--- Randomization test (Day 36)----
#' @title Prueba de Aleatorización y Visualización para el Análisis de Paralelismo Genómico (Día 36)
#'
#' @description Este script realiza una prueba de aleatorización para evaluar la significancia de la superposición en mutaciones adquiridas entre poblaciones, basada en los coeficientes de similitud de Sørensen-Dice. Construye un estadístico de prueba, realiza aleatorización y re-muestreo, calcula el valor p, y visualiza los resultados usando un diagrama de caja.
#'
#' @details 
#' El script realiza los siguientes pasos:
#' \itemize{
#'   \item Construye un estadístico de prueba que representa la diferencia en los coeficientes de similitud de Sørensen-Dice entre pares de muestras con el mismo gen de resistencia y aquellos con genes de resistencia diferentes para el Día 36 (P12).
#'   \item Realiza una prueba de aleatorización permutando las etiquetas ("mismo" o "diferente" gen de resistencia) y recalcula el estadístico de prueba 10,000 veces para generar una distribución del estadístico de prueba bajo la hipótesis nula.
#'   \item Calcula el valor p como la proporción de estadísticos de prueba permutados que son más extremos que el estadístico de prueba observado.
#'   \item Prepara los datos para la visualización y crea un diagrama de caja para ilustrar la superposición en mutaciones adquiridas entre muestras con el mismo o diferente gen de resistencia.
#' }
#'
#' @param sim_coef_NS_noP0 Un data frame que contiene coeficientes de similitud de Sørensen-Dice para cada par de muestras, anotado con información de pasaje y gen de resistencia.
#' @return Un objeto ggplot que muestra los coeficientes de similitud de Sørensen-Dice para pares de muestras con el mismo o diferente gen de resistencia, sin anotación del valor p.
#' @examples
#' \dontrun{
#' # Ejecutar el script para realizar la prueba de aleatorización y generar el gráfico
#' randomization_test_day36()
#' }
#' @import ggplot2
#' @importFrom stats sd
#' @export
#' 

# Construcción del Estadístico de Prueba
tmp2 <- sim_coef_NS_noP0[sim_coef_NS_noP0$passage1 == "P12" & 
                           sim_coef_NS_noP0$passage2 == "P12" & 
                           !is.na(sim_coef_NS_noP0$gene1) & 
                           !is.na(sim_coef_NS_noP0$gene2), ]

test_same_P12 <- tmp2[tmp2$same_gene == "same", "sim_coef"]
test_diff_P12 <- tmp2[tmp2$same_gene == "different", "sim_coef"]
test_gene_P12 <- mean(test_same_P12) - mean(test_diff_P12)

# Aleatorización y Re-muestreo
perm_gene_P12 <- c()
set.seed(123)
for (i in seq(1, 10000)) {
  tmp2$same_gene <- sample(tmp2$same_gene, replace = FALSE)
  perm_same_P12 <- mean(tmp2[tmp2$same_gene == "same", "sim_coef"])
  perm_diff_P12 <- mean(tmp2[tmp2$same_gene == "different", "sim_coef"])
  perm_gene_P12 <- c(perm_gene_P12, perm_same_P12 - perm_diff_P12)
}

# Valor p (proporción de permutaciones más extremas que el valor observado)
p_value_P12 <- length(perm_gene_P12[perm_gene_P12 >= test_gene_P12]) / length(perm_gene_P12)
print(p_value_P12)

# Preparación de Datos para la Visualización
test_gene_df <- data.frame(
  sim_coef = c(test_same_P12, test_diff_P12),
  label = c(rep("Same\nresistance\ngene", length(test_same_P12)), rep("Different\nresistance\ngene", length(test_diff_P12)))
)
colnames(test_gene_df) <- c("sim_coef", "label")
test_gene_df$label <- factor(test_gene_df$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene"))

#--- Fig 5A ---- 
# Creación y visualización del gráfico sin anotación del valor p
p <- ggplot(test_gene_df, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Overlap in acquired mutations\n(Paired Sørensen-Dice similarity)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55)

print(p)


## Identificación de genes mutados en diferentes poblaciones
#' @title Identificación y Análisis de Mutaciones de Genes en Diferentes Poblaciones
#'
#' @description Este script identifica qué genes están mutados en diferentes poblaciones y analiza las proporciones de mutaciones para cada gen dentro de varios grupos de genes de resistencia. También construye un gráfico para visualizar el porcentaje de poblaciones con mutaciones para cada gen.
#'
#' @details
#' El script realiza los siguientes pasos:
#' \itemize{
#'   \item Construye un data frame con los principales genes mutados en el pasaje P12.
#'   \item Crea una tabla de análisis de genes para el pasaje P12.
#'   \item Define vectores para genes de resistencia y anota cada población en función de estos vectores.
#'   \item Calcula el porcentaje de poblaciones con cada gen de resistencia que mutó en cada uno de los otros genes.
#'   \item Ordena y procesa el data frame para la visualización.
#'   \item Genera un gráfico de mapa de calor (Fig 5B) para visualizar las proporciones de mutación en diferentes genes principales y grupos de genes de resistencia.
#' }
#'
#' @return Un objeto ggplot que muestra el porcentaje de poblaciones con mutaciones para cada gen en diferentes grupos de genes de resistencia.
#'
#' @import ggplot2
#' @import viridis
#' @export
#'
#' @examples
#' # Cargar las bibliotecas necesarias
#' library(ggplot2)
#' library(viridis)
#'
#' # Suponiendo que 'breseq_NS_noP0' está cargado en el entorno

top_genes_in_study2 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(top_genes_in_study2) <- c("description", "Line", "hits")
top_genes <- names(tail(sort(table(top_genes_in_study2[top_genes_in_study2$hits > 0, "description"])), 20))

### Construcción de una tabla de análisis de genes para P12
breseq_gene_analysis2 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(breseq_gene_analysis2) <- c("gene", "Line", "hits")

# Definir los vectores de genes de resistencia
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  

# Anotar cada población según su gen de resistencia
breseq_gene_analysis2[breseq_gene_analysis2$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis2[breseq_gene_analysis2$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis2[breseq_gene_analysis2$Line %in% glycoside, "res_gene"] <- "glycoside"

# Eliminar filas donde res_gene es NA
breseq_gene_analysis2 <- breseq_gene_analysis2[!is.na(breseq_gene_analysis2$res_gene), ]

### Calcular el porcentaje de poblaciones con cada gen de resistencia que mutó en cada uno de los otros genes
breseq_gene_props <- data.frame(matrix(nrow = 0, ncol = 3))
l <- list(rfbA = rfbA, glycosyl = glycosyl, glycoside = glycoside)
names(l) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants")

# Bucle a través de genes y calcular proporciones de mutación
for (gene in top_genes) {
  for (j in seq_along(l)) {
    res_gene <- l[[j]]
    prop <- nrow(breseq_gene_analysis2[breseq_gene_analysis2$Line %in% res_gene & 
                                         breseq_gene_analysis2$gene == gene & 
                                         breseq_gene_analysis2$hits > 0, ]) / 
      nrow(breseq_gene_analysis2[breseq_gene_analysis2$Line %in% res_gene & 
                                   breseq_gene_analysis2$gene == gene, ])
    breseq_gene_props <- rbind(breseq_gene_props, c(gene, names(l)[j], prop))
  }
}

# Ordenar la columna 'res_gene' para que "rfbA mutants" aparezca primero, seguido por "PSPTO_4988 mutants" y luego "PSPTO_4991 mutants"
colnames(breseq_gene_props)
colnames(breseq_gene_props) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props$res_gene <- factor(breseq_gene_props$res_gene, levels = c("PSPTO_4991 mutants", "PSPTO_4988 mutants", "rfbA mutants"))

# Asignar nombres de columnas
colnames(breseq_gene_props) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props$percent_pops <- as.numeric(breseq_gene_props$percent_pops)
breseq_gene_props[breseq_gene_props$percent_pops == 0, "percent_pops"] <- 10^-4

# Reemplazar guiones problemáticos y asegurar codificación adecuada
breseq_gene_props$top_gene <- gsub("‑", "-", breseq_gene_props$top_gene)

# Mantener el orden de factores en el eje X
breseq_gene_props$top_gene <- factor(breseq_gene_props$top_gene, levels = rev(unique(breseq_gene_props$top_gene)))

#--- Plot Fig 5B----
ggplot(breseq_gene_props, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")) +
  scale_fill_viridis(name = "Populations\nwith mutation (%)", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  ylab("") +
  xlab("") +
  ggtitle("Day 36 of experimental evolution")



sessionInfo()
#R version 4.4.0 (2024-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 22631)

#Matrix products: default


#locale:
# [1] LC_COLLATE=Spanish_Mexico.utf8  LC_CTYPE=Spanish_Mexico.utf8    LC_MONETARY=Spanish_Mexico.utf8
#[4] LC_NUMERIC=C                    LC_TIME=Spanish_Mexico.utf8    

#time zone: America/Mexico_City
#tzcode source: internal

#attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] viridis_0.6.5     viridisLite_0.4.2 ggplot2_3.5.1     readxl_1.4.3     

#loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.3         knitr_1.48        rlang_1.1.4       xfun_0.47        
#[6] generics_0.1.3    labeling_0.4.3    glue_1.7.0        colorspace_2.1-0  htmltools_0.5.8.1
#[11] gridExtra_2.3     scales_1.3.0      fansi_1.0.6       rmarkdown_2.28    grid_4.4.0       
#[16] cellranger_1.1.0  evaluate_0.24.0   munsell_0.5.1     tibble_3.2.1      fastmap_1.2.0    
#[21] lifecycle_1.0.4   compiler_4.4.0    dplyr_1.1.4       pkgconfig_2.0.3   rstudioapi_0.16.0
#[26] farver_2.1.2      digest_0.6.37     R6_2.5.1          tidyselect_1.2.1  utf8_1.2.4       
#[31] pillar_1.9.0      magrittr_2.0.3    withr_3.0.1       tools_4.4.0       gtable_0.3.5     
