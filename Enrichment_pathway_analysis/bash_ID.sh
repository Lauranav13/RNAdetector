#Código en bash para sacar los IDs coincidentes tanto de los genes totales, como de los genes diferencialmente expresados. 

cat ID_TABLE.txt | tr ',' '\n' > id_table.txt


#Ordenamos y quitamos duplicados de cada archivo.
sort -u id_table.txt -o id_table_sorted.txt
sort -u all_genes.txt -o all_genes_sorted.txt

#Comparamos los dos archivos.
common=$(comm -12 id_table_sorted.txt all_genes_sorted.txt)

#Imprimimos por pantalla los valores coincidentes.
echo "Common values:"
echo "$common"

#Imprimimos por pantalla la cantidad de valores que coinciden.
common_count=$(echo "$comunes" | wc -l)
echo "Cantidad de valores comunes: $common_count"
#Ordenamos y quitamos duplicados.
sort -u differentially_expressed_genes.txt -o differentially_expressed_genes_sorted.txt

#Comparamos los dos archivos.
common=$(comm -12 id_table_sorted.txt differentially_expressed_genes_sorted.txt)


#Imprimimos por pantalla los valores coincidentes.
echo "Common values:"
echo "$common"

#Imprimimos por pantalla la cantidad de valores que coinciden.
common_count=$(echo "$common" | wc -l)
echo "Cantidad de valores comunes: $common_count"


#Después de realizar esto, podemos correr el script de python pathway.py y calcular el valor de la hipergeométrica en cada una de las rutas con el script de perl hypergeom.pl.
