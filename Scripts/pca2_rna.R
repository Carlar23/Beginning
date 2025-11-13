###########################################
# Análisis de Componentes Principales (PCA)
###########################################

# Directorio de trabajo
setwd("C:/Users/diazc/Desktop/Algorit_inteligen_artificial/tema_1/data_export/rna_cancer")

# Librerías necesarias
library(ggplot2)
library(dplyr)

# Lectura de los datos
datos_raw <- read.csv("data.csv")
etiquetas <- read.csv("labels.csv")

# Seleccionamos los primeros 500 genes (ignorando la primera columna si es ID)
genes <- datos_raw[, 2:501] %>% mutate_all(as.numeric)

# Inspección básica de datos
cat("Dimensiones del dataset:", dim(genes), "\n")
cat("Valores NA presentes:", anyNA(genes), "\n")

# Conteo de ceros y NA
resumen_calidad <- data.frame(
  Variable = names(genes),
  NAs = colSums(is.na(genes)),
  Ceros = colSums(genes == 0)
)

# Mostrar las 10 variables con más ceros
head(resumen_calidad[order(-resumen_calidad$Ceros), ], 10)

# Gráfico de ceros por variable
ggplot(resumen_calidad, aes(x = reorder(Variable, -Ceros), y = Ceros)) +
  geom_col(fill = "steelblue") +
  labs(title = "Cantidad de ceros por variable", x = "Variable", y = "Cantidad de ceros") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "none")

# Boxplot de los primeros 10 genes
boxplot(genes[, 1:10], main = "Distribución de los primeros 10 genes", col = "lightblue")

# ---------------------------------------------------
# Cálculo del PCA
# ---------------------------------------------------

pca_modelo <- prcomp(genes, center = TRUE, scale. = TRUE)  # Aquí escalamos las variables

# Varianza explicada
var_explicada <- (pca_modelo$sdev)^2 / sum(pca_modelo$sdev^2)
var_acumulada <- cumsum(var_explicada)

# Número mínimo de componentes que explican el 90% de la varianza
num_pc <- which(var_acumulada >= 0.90)[1]
cat("Número de componentes que explican al menos el 90% de la varianza:", num_pc, "\n")

# DataFrame con los resultados del PCA
pca_df <- as.data.frame(pca_modelo$x)
pca_df$Clase <- etiquetas$Class

# Etiquetas de ejes dinámicas
x_lab <- paste0("PC1 (", round(var_explicada[1] * 100, 1), "%)")
y_lab <- paste0("PC2 (", round(var_explicada[2] * 100, 1), "%)")

# ---------------------------------------------------
# Visualización de los resultados del PCA
# ---------------------------------------------------

grafico_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Clase)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA - Tipos de Cáncer", x = x_lab, y = y_lab, color = "Clase") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# Mostrar el gráfico
print(grafico_pca)

# ---------------------------------------------------
# Guardar el gráfico
# ---------------------------------------------------

ggsave("grafico_pca_v2.pdf", plot = grafico_pca, width = 8, height = 6)

# ---------------------------------------------------
# Gráfico adicional: varianza explicada
# ---------------------------------------------------

df_varianza <- data.frame(
  Componente = paste0("PC", 1:length(var_explicada)),
  Varianza = var_explicada,
  VarianzaAcumulada = var_acumulada
)

ggplot(df_varianza, aes(x = Componente, y = Varianza)) +
  geom_col(fill = "darkorange") +
  geom_line(aes(y = VarianzaAcumulada), group = 1, color = "black", size = 1) +
  geom_point(aes(y = VarianzaAcumulada), color = "black", size = 2) +
  labs(title = "Varianza explicada por cada componente principal",
       x = "Componentes principales", y = "Proporción de varianza") +
  theme_minimal()

