###############################
# Grafico varianza n-dim MDS  #
###############################

# Directorio en donde tendremos nuestros datos
setwd("~/Desktop/PER_13711/Scripts/data/rna_cancer")

# Calcular la matriz de distancias en tus datos
datos <- as.matrix(read.csv("datos_500.csv", row.names = 1))
distancia <- dist(datos)

# Aplicar MDS a un número alto de dimensiones (por ejemplo, 10)
mds_result <- cmdscale(distancia, k = 10, eig = TRUE)

# Calcular la varianza explicada
# La varianza explicada para cada dimensión se calcula usando los valores propios (eigenvalues)
valores_propios <- mds_result$eig
varianza_explicada <- valores_propios / sum(valores_propios)
varianza_acumulada <- cumsum(varianza_explicada)

# Graficar la varianza explicada y la varianza acumulada
library(ggplot2)

# Crear data frame para la gráfica
df_varianza <- data.frame(
  Dimensión = 1:10,
  Varianza_Explicada = varianza_explicada[1:10],
  Varianza_Acumulada = varianza_acumulada[1:10]
)

# Gráfico de varianza explicada
ggplot(df_varianza, aes(x = Dimensión)) +
  geom_bar(aes(y = Varianza_Explicada), stat = "identity", fill = "skyblue") +
  geom_line(aes(y = Varianza_Acumulada), color = "red", size = 1) +
  geom_point(aes(y = Varianza_Acumulada), color = "red") +
  labs(title = "Varianza Explicada y Acumulada en MDS",
       x = "Dimensión",
       y = "Varianza Explicada") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Varianza Acumulada")) +
  theme(axis.title.y.right = element_text(color = "red"))



###############################
# Multidimensional Scaling 3D #
###############################


# Instalar plotly si no lo tienes instalado
# install.packages("plotly")

library(plotly)

# Supongamos que 'datos' es tu matriz de expresión génica
# Asegúrate de que los datos están en formato numérico y sin columnas de etiquetas

# Calcular la matriz de distancias
distancia <- dist(datos)

# Aplicar MDS a 3 dimensiones
mds_result <- cmdscale(distancia, k = 3)

# Convertir el resultado a un data frame
mds_df <- as.data.frame(mds_result)
colnames(mds_df) <- c("Dim1", "Dim2", "Dim3")

# Añadir etiquetas de cáncer
labels <- read.csv("labels.csv")  # Asegúrate de cambiar "ruta" por la ubicación de tu archivo
mds_df$Cancer_Type <- labels$Class

# Graficar en 3D
fig <- plot_ly(mds_df, x = ~Dim1, y = ~Dim2, z = ~Dim3, 
               color = ~Cancer_Type, colors = "Set1", 
               marker = list(size = 4))

# Mostrar el gráfico
fig
