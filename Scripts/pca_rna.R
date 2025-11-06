#######################################
# Analisis de componentes principales #
#######################################

# Directorio en donde tendremos nuestros datos
setwd("~/Desktop/PER_13711/Scripts/data/rna_cancer")

# Carga de librerias
library(stats)   # librería para el PCA
library(ggplot2) # librería para hacer la representación gráfica

# Lectura de datos
data.raw <- read.csv('data.csv')
labels.raw <- read.csv('labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- data.frame(sapply(data.raw[2:501], as.numeric))

# Estadisticos
summary(data[, 1:10])

# Valores NA y 0
anyNA(data)
na_counts <- colSums((is.na(data)))

any(data == 0)
zero_counts <- colSums(data == 0)

# Graficamos
zero_df <- data.frame(
  Variable = names(zero_counts),
  Zeros = as.numeric(zero_counts)
)
ggplot(zero_df, aes(x = Variable, y = Zeros, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Cantidad de ceros por columna",
       x = "Variable",
       y = "Número de ceros") +
  theme_minimal() +
  theme(legend.position = "none")  # Oculta la leyenda

# Podemos hacer un diagrama de cajas para variable y vemos los estadísticos y outliers
boxplot(data[, 1:10], main = "Boxplot de los 10 primeros genes")


# Funcion prcomp()
#   data: conjunto de datos
#   center: si queremos que las variables esten centradas en cero
#   scale: si queremos que las variables tengan varianza 1



# Calculo de componentes principales con la funcion prcomp
pca.results <- prcomp(data, center=TRUE, scale.=FALSE)

# Resultado de las componentes principales
pca.df <- data.frame(pca.results$x)

# Varianza (cuadrado de la desviacion tipica)
varianzas <- pca.results$sdev^2

# Total de la varianza de los datos
total.varianza <- sum(varianzas)

# Varianza explicada por cada componente principal
varianza.explicada <- varianzas/total.varianza

# Calculamos la varianza acumulada 
varianza.acumulada <- cumsum(varianza.explicada)

# Tomamos el numero de componentes principales que explican el 90% de la varianza
n.pc <- min(which(varianza.acumulada > 0.90))

# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')

# Representación gráfica de las primeras dos componentes principales respecto a los datos
ggplot(pca.df, aes(x=PC1, y=PC2, color=labels.raw$Class)) +
  geom_point(size=3) +
  scale_color_manual(values=c('red', 'blue', 'green', 'orange', 'purple')) +
  labs(title='PCA - Types of Cancer', x=x_label, y=y_label, color='Grupo') +
  theme_classic() +
  theme(panel.grid.major = element_line(color="gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5))
