# Analysis of bacterial genomes evolved in the absence of phage
# Análisis de genomas bacterianos evolucionados en ausencia de fagos
# Reena Debray
# Feb 1, 2022
# Modified: Evelia Coss (03/Sep/2024)
# Paper: https://academic.oup.com/mbe/article/39/9/msac182/6673247?login=false#371728839
# Github: https://github.com/reenadebray/loss-of-resistance/tree/main
# Script: https://github.com/reenadebray/loss-of-resistance/blob/main/Phenotype_data_analysis.R

# --- Paquetes ----
library(growthrates) # fit_growthmodel(): Calcular la tasa de crecimiento a partir de datos de microplacas.
# NOTA: Toma los datos de la siguiente forma: cada fila es la medida de un único pocillo en un único punto temporal, 
# con columnas adicionales «Tratamiento» y «DO». Se ajusta una curva logística a cada pocillo por separado utilizando 
# la DO a lo largo del tiempo. 
# He utilizado los siguientes parámetros: lower=c(y0=0.000001,mumax=0,K=0),upper=c(y0=0.05,mumax=5,K=1.5)),p=c(y0=0.01,mumax=0.2,K=0.1)
library(lmerTest) # Pruebas estadisticas
# Paquetes adicionales (no los mencionan pero son necesarios)
library(readxl) # lectura de archivos de excel
library(tidyverse) # Manipulacion de datos
library(reshape2) # melt(): reduccion de datos

#---- GC_long_to_growthrates() -------

#' @title Calcular tasas de crecimiento
#' @description
#' Esta función procesa los datos de crecimiento de un experimento en múltiples pozos, ajusta un modelo logístico de 
#' crecimiento para cada pozo, y extrae la tasa de crecimiento r. 
#' @author Reena Debray
#' 
#' @param GC_long # Un data frame en formato largo que contiene los datos de crecimiento.
#' @param lower   # Un vector de límites inferiores para los parámetros del modelo de crecimiento.
#' @param upper   # Un vector de límites superiores para los parámetros del modelo de crecimiento.
#' @param p       # Un vector de parámetros iniciales para el ajuste del modelo.
#'
#' @return # El resultado es un data frame con las tasas de crecimiento para cada combinación de pozo y tratamiento.
#' \item{well}{The identifier of the well.}
#' \item{treatment}{The treatment associated with the well.}
#' \item{r}{The growth rate obtained from the fitted logistic model.}
#' 
#' @import growthrates 
#' @import lmerTest
#'
#' @export
#' 
#' @examples
#' # Example data
#' GC_long <- data.frame(
#'   Well = rep(c("A1", "A2"), each = 10),
#'   Time = rep(1:10, 2),
#'   OD = runif(20),
#'   Treatment = rep(c("Control", "Treatment"), each = 10)
#' )
#'
#' # Parameters for the logistic model
#' lower <- c(y0 = 0.000001, mumax = 0, K = 0)
#' upper <- c(y0 = 0.05, mumax = 5, K = 1.5)
#' p <- c(y0 = 0.01, mumax = 0.2, K = 0.1)
#'
#' # Calculate growth rates
#' growth_rates <- GC_long_to_growthrates(GC_long, lower, upper, p)
#' print(growth_rates)

GC_long_to_growthrates <- function(GC_long, lower, upper, p){
  ### initialize data frame / Iniciar dataframe
  # Objeto donde se almacenarán las tasas de crecimiento calculadas
  growthrates<-data.frame(matrix(nrow=0, ncol=3)) 
  # Populate with model fit / Rellenar con el ajuste del modelo
  # Argumentos:
  # `treatment`: Para cada pozo, se extrae el tratamiento asociado a partir de la columna Treatment.
  # Modelo de crecimiento: Se ajusta un modelo de crecimiento logístico (especificado por `grow_logistic`) a los datos de tiempo (`Time`) 
  # y densidad óptica (`OD`) para el pozo actual.
  # `r`: Se extrae el coeficiente de crecimiento (r) del modelo ajustado usando coef, y se convierte a un número con as.numeric.
  # Agregar a growthrates: Se añade una nueva fila al data frame growthrates que contiene el pozo (well), 
  # el tratamiento (treatment) y la tasa de crecimiento (r).
  for (well in unique(GC_long$Well)){
    treatment <- GC_long[GC_long$Well == well,"Treatment"][1]
    r <- as.numeric(coef(fit_growthmodel(FUN = grow_logistic, GC_long[GC_long$Well==well,"Time"], GC_long[GC_long$Well==well,"OD"],p=p,lower=lower,upper=upper))[2])
    # Cambiar formato a numerico
    growthrates <- rbind(growthrates,c(well,treatment,r))
  }
  ### Return output
  colnames(growthrates)=c("well","treatment","r")
  growthrates$r <- as.numeric(growthrates$r)
  return(growthrates)
}

## ---- Figure 1: Resistance is costly -------

### Read in "Costs_of_Res.xlsx" as costs_of_res / Cargar informacion de "Costs_of_Res.xlsx" como variable costs_of_res
costs_of_res <- read_excel("Costs_of_Res.xlsx")
### Calculate fitted values controlling for phage resistance and plate layout / Calcular los valores ajustados controlando la resistencia del fago y la disposición de la placa
### Express fitted values as a percentage of wild-type fitness / Expresar los valores ajustados como porcentaje de la aptitud de tipo wild-type
costs_of_res$fitted <- fitted.values(lm(r~(Population=="ancDC3000") + Column + Plate,costs_of_res))
costs_of_res$fitted_percWT <- costs_of_res$fitted/mean(unlist(costs_of_res[costs_of_res$Population=="ancDC3000","fitted"]))*100
### Reorder the isolates by resistance gene and growth rate / Reordenar los aislados por gen de resistencia y tasa de crecimiento
costs_of_res<-costs_of_res[order(costs_of_res$Gene,-costs_of_res$r),]
# place the isolate with no detected genetic differences on the right-hand side of the graph / Coloque el aislado sin diferencias genéticas detectadas en la parte derecha del gráfico
costs_of_res[costs_of_res$Population=="MS15","Gene"] <- "Z" 
# place the isolate with no detected genetic differences on the right-hand side of the graph / Coloque el aislado sin diferencias genéticas detectadas en la parte derecha del gráfico
costs_of_res[costs_of_res$Population=="MS15","Annotation"] <- "Z" 
costs_of_res$order<-seq(1,nrow(costs_of_res))

### Resistant bacteria grow more slowly than their sensitive ancestor / Las bacterias resistentes crecen más despacio que sus antepasadas sensibles
costs_of_res_R<-costs_of_res[costs_of_res$Population!="ancDC3000",]
costs_agg<-aggregate(costs_of_res_R$fitted_percWT,by=list(costs_of_res_R$Population,costs_of_res_R$Gene),FUN=mean)
colnames(costs_agg)<-c("Population","Gene","fitted_percWT")
t.test(x = costs_agg$fitted_percWT,mu = 100,alternative = "less")
### Variation in growth rates is not explained by resistance gene (exclude population with no detected mutations) / La variación en las tasas de crecimiento no se explica por el gen de resistencia (excluir la población sin mutaciones detectadas)
anova(lm(fitted_percWT~Gene,costs_agg[costs_agg$Population!="MS15",]))

### ---- Figura 1B. Genetic mechanisms and fitness costs of phage resistance / Mecanismos genéticos y costes de fitness de la resistencia a los fagos -----
# Population growth rates of resistant strains in the absence of phage, based on a logistic model fitted to a 40 h growth curve (n = 22 strains)
# Tasas de crecimiento de la población de cepas resistentes en ausencia de fago, basadas en un modelo logístico ajustado a una curva de crecimiento de 40 h (n = 22 cepas).
ggplot(costs_of_res[costs_of_res$Population!="ancDC3000",]) + 
  stat_summary(aes(reorder(Population,-fitted_percWT),fitted_percWT,group=Population,fill=Gene),geom="pointrange",shape=21,size=0.8) + 
  theme_classic(base_size=18) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(fill=F) + 
  facet_grid(~paste(Gene,Annotation,sep="\n"), scales="free_x",space="free_x") + 
  theme(strip.background = element_blank(),strip.placement="outside") + 
  scale_fill_brewer(palette="Dark2") + 
  ylab("Fitness (% of wild-type)") + 
  geom_hline(yintercept=100,linetype="dashed")

## ---- Figure 2: Phenotypic changes in bacterial populations during experimental evolution in the absence of phages / Cambios fenotípicos en poblaciones bacterianas durante la evolución experimental en ausencia de fagos. -------

### Read in "Fitness_over_Time.xlsx" as fitness_over_time / Cargar informacion de "Fitness_over_Time.xlsx" como variable fitness_over_time
fitness_over_time <- read_excel("Fitness_over_Time.xlsx")
### Calculate fitted values controlling for phage resistance, passage, and plate layout / Calcular los valores ajustados controlando la resistencia del fago, el paso y la disposición de la placa.
fitness_over_time$fitted <- fitted.values(lm(r~Passage*Type+Column,fitness_over_time))
### Express fitted values as a percentage of wild-type fitness / Expresar los valores ajustados como porcentaje de la desarrollo de tipo wild-type
for (passage in seq(0,12,2)){
  anc_means<-mean(unlist(fitness_over_time[fitness_over_time$Passage==passage & fitness_over_time$Type=="ANC","fitted"]))
  fitness_over_time[fitness_over_time$Passage==passage & fitness_over_time$Type=="PR","fitted_percWT"]<-100*fitness_over_time[fitness_over_time$Passage==passage & fitness_over_time$Type=="PR","fitted"]/anc_means
}

### ----- Figure 2A: Growth over time / Crecimiento a lo largo del tiempo -----
df <- data.frame(aggregate(unlist(fitness_over_time[fitness_over_time$Type=="PR","fitted_percWT"]), 
                            by=list(unlist(fitness_over_time[fitness_over_time$Type=="PR","Passage"])),FUN=mean), 
                 aggregate(unlist(fitness_over_time[fitness_over_time$Type=="PR","fitted_percWT"]),
                            by=list(unlist(fitness_over_time[fitness_over_time$Type=="PR","Passage"])),FUN=sd)[,2])
colnames(df) <- c("passage","mean","SD")

# Grafica
ggplot() + geom_point(data=df,aes(x=passage*3,y=mean),size=3) +
  geom_errorbar(data=df,aes(x=passage*3,ymin=mean-SD,ymax=mean+SD),width=0,size=1) + 
  geom_line(data=df,aes(x=passage*3,y=mean)) + 
  theme_classic(base_size=16) + 
  geom_hline(yintercept=100,linetype="dashed") + 
  scale_x_continuous(breaks=seq(0,36,6)) + 
  ylab("Fitness (% of wild-type)") + xlab("Day of experimental evolution") + 
  coord_cartesian(ylim=c(60,115))


### ---- Figure 2B: Resistance over time / Resistencia en el tiempo ------

# NOTA: Abri el excel y corre un error de compatibilidad.
### Read in "Resistance_over_Time.xlsx" as res_over_time / Cargar informacion de "Resistance_over_Time.xlsx" como variable res_over_time
res_over_time <- read_excel("Resistance_over_Time.xlsx", sheet = "Data")
# Transformación de los datos en formato largo
res_long <- melt(res_over_time, id.vars = "isolate", measure.vars = c("Prop_S_P0","Prop_S_P4","Prop_S_P8","Prop_S_P12"))
res_long$variable <- as.character(res_long$variable) # Conversión de la columna variable a caracteres:  
res_long$passage <- as.numeric(substr(res_long$variable,9,nchar(res_long$variable)))
res_long <- res_long[!res_long$isolate %in% c("FMS6","SNK6"),]  # Filtrado de aislamientos


### ----- Code for Figure 2B  ---------

### Excludes populations with colonies that appeared sensitive on plates but were not affected by phage in growth curves
# Excluye las poblaciones con colonias que parecían sensibles en las placas pero que no se veían afectadas por el fago en las curvas de crecimiento

ggplot(res_long, aes(passage*3,value*100,group=isolate,color=(isolate%in%c("SNK7","QAC5","VCM4")))) +
  geom_point(size=4.5,shape=18) +
  geom_line(linewidth=0.8)+
  theme_classic() +
  guides(color=F)+
  scale_color_brewer(palette="Dark2") + 
  scale_x_continuous(breaks=c(0,12,24,36)) + 
  theme_classic(base_size=16) + 
  ylab("Proportion of phage sensitivity in population (%)") +
  xlab("Day of experimental evolution") + 
  theme(axis.title.y=element_text(size=14))

### Increase in fitness over time; no main or interaction effect of populations that reverted to sensitivity / 
# Aumento de la aptitud con el tiempo; sin efecto principal o de interacción de las poblaciones que revirtieron a la sensibilidad.
fitness_over_time[fitness_over_time$Population %in% c("SNK7","QAC5","VCM4"),"Outcome"] <- "S"
fitness_over_time[!fitness_over_time$Population %in% c("SNK7","QAC5","VCM4"),"Outcome"] <- "R"
anova(lmer(fitted_percWT~Passage*Outcome+(1|Population),fitness_over_time))

### No effect of costs of resistance on reversion to sensitivity / Sin efecto de los costes de resistencia sobre la reversión a la sensibilidad
costs_agg[costs_agg$Population%in%c("SNK7","QAC5","VCM4"),"Outcome"] <- "S"
costs_agg[!costs_agg$Population%in%c("SNK7","QAC5","VCM4"),"Outcome"] <- "R"
t.test(fitted_percWT~Outcome, costs_agg)


## ----- Figure 3: Replay experiment / Experimento de repetición ---------
### Read in "Replay.xlsx" as replay / Cargar informacion de "Replay.xlsx" como variable replay
replay <- read_excel("Replay.xlsx", sheet = "Data")
replay$Prop_S_P12 <- as.numeric(as.character(replay$Prop_S_P12)) # Convert as numeric / Convertir a numerico

# NOTA: no encuentro Prop_R_P12, es mas bien Prop_S_P12
### Populations derived from a founder that re-evolved sensitivity were more likely to re-evolve sensitivity as well
# Las poblaciones derivadas de un fundador que reevolucionó la sensibilidad tenían más probabilidades de reevolucionar también la sensibilidad
t.test(Prop_S_P12~Founder_outcome, data=replay)

### ----- Figure 3C: Phage resistance outcomes in the replay experiment mirror the outcomes of their founding population ---------
# Proportion of colonies in each population that were scored as phage-sensitive after 36 days in the absence of phages, 
# as indicated by disruption in bacterial growth upon encountering phage on an agar plate (n=54 populations with 96 colonies 
# sampled per population).
agg <- data.frame(
  aggregate(100*replay$Prop_S_P12, by=list(replay$Founder_outcome), FUN=function(x){mean(x,na.rm=T)}) ,
  aggregate(100*replay$Prop_S_P12,by=list(replay$Founder_outcome),FUN=function(x){sd(x,na.rm=T)/sqrt(length(x))}))

# Rename
colnames(agg) <- c("Founder_outcome","mean","Founder_2","SE")
ggplot() +
  geom_errorbar(data=agg,aes(Founder_outcome,ymax=mean+SE,ymin=mean-SE),alpha=0.4,width=0.6,size=0.8) +
  stat_summary(data=agg,aes(Founder_outcome,mean,fill=Founder_outcome),geom="bar",fun=mean,width=0.8,color="black") +
  geom_jitter(data=replay,aes(Founder_outcome,Prop_S_P12*100),size=3,alpha=0.4,width=0.05,height=1) + 
  theme_classic(base_size=16) +
  scale_fill_brewer(palette="Dark2") + 
  guides(fill=F) + 
  xlab("") + ylab("Proportion of phage sensitivity in replay populations (%)") + 
  theme(axis.text.x=element_blank())


### ---- Supplementary figure / Figura suplementaria ------

# diminishing returns / rendimientos decrecientes
cast <- dcast(fitness_over_time[fitness_over_time$Sample!="ANC1_P0",],Population~Passage,FUN=mean,value.var = "fitted_percWT")
colnames(cast) <- c("Population","P0","P2","P4","P6","P8","P10","P12")
# Analisis de ANOVA
anova(lm((P12-P0)~P0,cast))

# Graph / Grafica
ggplot(cast,aes(P0,P12-P0)) + 
  stat_smooth(method="lm",color="gray30") + 
  geom_jitter(width=0.8,height=0.8,size=2,alpha=0.8) + 
  theme_classic(base_size=16) + 
  xlab("Initial fitness\n(% wild-type, day 0)") + ylab("Change in fitness\n(% wild-type, day 0 - day 36)")


