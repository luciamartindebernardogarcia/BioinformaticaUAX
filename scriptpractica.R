#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 22 OCTUBRE 23:59
## Se requiere la entrega de este script completado con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #esto mide las dimensiones
head(data) #así mostramos las primeras 
tail(data) #muestra las últimas filas
# Hacemos un primer histograma para explorar los datos
hist(data)#si lo queremos guardar, simplemente clicamos la foto
hist (data, col = "pink", main="GSE5583 - Histogram")
# Transformamos los datos con un logaritmo 
data_log=log2(data)
hist(data_log)
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
Sirve para poner de forma más adecuada los datos y que así queden mejor agrupados y por lo tanto con una distribución normal.

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
boxplot(data_log)
boxplot(data_log,col=c("blue","blue","blue","orange","orange","orange"))
boxplot(data_log,col=c("blue","blue","blue","orange","orange","orange"), main="GSE5583-boxplots",las=2)
#las sirve para poner los datos en vertical, col para poner colores y main para nombrar el boxplot
# ¿Qué es un boxplot?
Un diagrama de barras o de bigotes y nos representa la mediana y los cuartiles

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
hc = hclust(as.dist(1-cor(data_log)))
hc es el clustering y después hemos usado plot para representarlo
plot(hc, main="GSE5583 - Hierarchical Clustering")
Sí, la separación es correcta
#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt <- data[,1:3]
#la coma delante separa las filas (1) de las columnas (3) 
ko <- data [,4:6]
class (wt)
Esto sirve para separar las columnas de una tabla en WT y KO
Los datos que hemos generado son una tabla o matriz
head (wt)#sirve para ver encabezado de las wt

# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply (wt, 1, mean)
Calcula la función que yo quiera con los datos que yo quiera
Si ponemos wt.mean nos sale cada gen y debajo su media
head (wt.mean)
ko.mean = apply (ko, 1, mean) #El 1 indica que queremos calcular la función para cada fila y el 2 para cada columna
head (ko.mean)
# ¿Cuál es la media más alta?
max(wt.mean)
max (ko.mean)
La media más alta corresponde a la de los KO
# Ahora hacemos un scatter plot (gráfico de dispersión)
plot (ko.mean ~ wt.mean) #El símbolo de aproximado indica el enfrentamiento entre los dos ejes
plot (ko.mean ~ wt.mean, xlab = "WT", ylab = "KO", main = "GSE5583 - Gráfico de dispersión")

# Añadir una línea diagonal con abline (queremos poner una línea de regresión. Se hace cuando el plot está abierto)
abline (0, 1, col = "red")
Teniendo en cuenta que y = mx+n, nuestra recta es igual que y = 1*x = x. 0 es n y 1 la m
Si queremos tener una recta horizontal ponemos: abline (h=2, col="blue")
Si queremos tener una recta vertical ponemos: abline (v=5,col="green")
# ¿Eres capaz de añadirle un grid?
grid()
Esto es para poner una cuadrícula
# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean
Obtenemos una tabla con las correspondientes medias
# Hacemos un histograma de las diferencias de medias
hist (diff.mean)
Obtenemos un histograma sin tener en cuenta los logaritmos por lo que sale más "feo", si lo queremos mejor podemos hacer:
Si queremos ponerle color al histograma: hist (diff.mean, col = "pink")
# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
Los tres apartados los tenemos a continuación: 
Por cada fila tenemos que hacer un t-test y por cada uno sacamos un p-value. 
pvalue = NULL
tstat = NULL
for (i in 1 : nrow (data)) { #Para cada gen
       x = wt[i,] #gene wt número i
       y = ko[i,] #gene ko número i

     # Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
 	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic

}
head(pvalue)
La i significa iterar: para cada fila desde la uno hasta el final de filas (nrow daba el nº de filas totales), haz lo que ponemos ahí. 
i solo la podemos llamar si está dentro del bucle
x = wt[i,] e y = ko[i,] son vectores. Solo ponemos i y no ponemos nada más porque queremos que coja todo
Hacemos el t-test y lo estans guardando en una variable que llamamos t

length (pvalue)

# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
Los usamos sin transofrmar porque estaríamos modificando nuestros resultados. A transformar podemos referirnos a hacer logaritmos 
# ¿Cuántas valores tiene cada condición?
En cada condición tenemos 3 valores. Cada gen tiene 6 muestras o réplicas (3 wt y 3 ko) y 2 condiciones (ko y wt)

# Ahora comprobamos que hemos hecho TODOS los cálculos
Lo comprobamos en casa. 

# Hacemos un histograma de los p-values.
hist (pvalue)
# ¿Qué pasa si le ponemos con una transformación de -log10?
La transformación la haremos porque nos interesan los pvalues pequeños (<0,05)
hist (-log10 (pvalue),col= "pink")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot (diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")
# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?

diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline (v = diff.mean_cutoff, col = "blue", lwd = 3)
abline (h = -log10 (pvalue_cutoff), col = "green", lwd = 3)

#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)Esto se pondría para el -2


Ponemos el log10 del pvalue porque al generar el volcano con logaritmos sino no nos lo reconocerá
lwd es el ancho de línea
# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
La diferencia de medias es 2 y -2. Abs es para hacer un valor absoluto
Para evitar tener que hacer un adicional para valores negativos, ponemos valor absoluto y sacamos todos aquellos por encima de 2
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff

# Ahora el filtro de p-value
Extraemos todos los pvalues que sean menores o iguales al filtro de pvalues que hemos puesto

filter_by_pvalue = pvalue <=pvalue_cutoff
dim(data[filter_by_pvalue, ])

Dim=dimensiones. Estoy buscando cuántos genes sobrepasan ese filtro. 

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
Combinamos los filtros y nos quedamos con los genes comunes a ambos filtros. Los genes que cumplen los dos criterios son 426 
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered =data[filter_combined,]
dim (filtered)
head (filtered)


# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo

plot (diff.mean, -log10 (pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]), col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
Cuando hemos calculado la diferencia de medias, hemos calculado los wt - ko. Los ko son los que están sobreexpresados por lo que si los valores de ko son mayores que los de wt, la resta sale con valor negativo. Por lo tanto, los reprimidos salen en el lado positivo y los sobreexpresados en el negativo. 
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0], 
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
heatmap (filtered)
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, labRow=FALSE)

Para las columnas hemos agrupado por wt o ko y las filas hemos agrupado como sobreexpresado o reprimido

# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
cexCol es el tamaño de la letra en el eje x
labRow=FALSE significa que hemos quitado los nombres de los genes
Colv y Rowv son los dendrogramas (indicamos filas y columnas)

# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
heatmap (filtered)

heatmap (filtered, Rowv=rowv, Colv=colv, cexCol=0.7, col=hcl.colors(50))


######A PARTIR DE AQUÍ NO ES NECESARIO REALIZARLO#######

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
La línea discontínua indica el 0. Veremos que los azules (reprimidos) estarán por debajo del cero y los rojos (sobreexpresados) por encima. 
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)


# Hacemos nuestro heatmap

heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")

# Lo guardamos en un archivo PDF

pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row", labRow=FALSE)
dev.off()
pdf ("GSE5583_DE_Heatmapverdeyrojo.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
col = redgreen(75), scale ="row",labRow=FALSE)
dev.off()

# Guardamos los genes diferencialmente expresados y filtrados en un fichero

write.table (filtered, "GSE5583_DE.txt", sep = "\t", 
	quote = FALSE)