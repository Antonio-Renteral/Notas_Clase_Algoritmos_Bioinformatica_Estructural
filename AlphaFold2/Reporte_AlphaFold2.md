Durante este ejercicio trabaje en el modelado de la estructura terciaria de una proteina utilizando AlphaFold2. Para ello
seleccione la secuencia identificada como B1AK53 - ESPN_HUMAN, correspondiente a la proteina Espin humana, la cual esta asociada
a la organizacion del citoesqueleto y estructuras celulares especializadas.

La secuencia de aminoacidos fue obtenida a partir de una base de datos publica y posteriormente utilizada como entrada en el
servidor de AlphaFold2 para predecir su estructura tridimensional. A partir de la secuencia primaria, el algoritmo genero un
modelo estructural completo sin necesidad de una plantilla experimental previa, basandose en patrones aprendidos de estructuras
conocidas.

Una vez obtenida la prediccion, se evaluo la calidad del modelo utilizando las metricas proporcionadas por AlphaFold2.
En particular, se analizo el valor de pLDDT como indicador de confianza local por residuo. Se observaron regiones con valores
altos de pLDDT, lo que sugiere una prediccion confiable en los dominios estructurados de la proteina. Tambien se identificaron
regiones con valores mas bajos, las cuales pueden corresponder a segmentos flexibles o potencialmente desordenados, algo comun en
proteinas con funciones regulatorias o de interaccion.

Ademas del pLDDT, se considero la coherencia general del plegamiento y la organizacion de posibles dominios estructurales dentro
del modelo generado.

Posteriormente, el modelo fue evaluado mediante SwissModel Assess, herramienta que permite analizar parametros de calidad
estructural como consistencia geometrica y compatibilidad entre secuencia y estructura. Los resultados indicaron que el modelo
presenta una calidad global adecuada y no muestra desviaciones estructurales significativas.

Como complemento, tambien se utilizo la plataforma SAVES para realizar una validacion adicional del modelo. Esta herramienta
integra distintos metodos de evaluacion que permiten examinar la calidad estereoquimica y detectar posibles regiones problematicas.
En general, los resultados obtenidos fueron consistentes con una estructura correctamente plegada en las regiones de mayor confianza.

En conclusion, este ejercicio permitio aplicar herramientas de prediccion estructural basadas en inteligencia artificial para
modelar la estructura terciaria de una proteina humana a partir de su secuencia de aminoacidos. Asimismo, el analisis de las
metricas de calidad fue fundamental para interpretar la confiabilidad del modelo y comprender que, aunque estas herramientas son
altamente precisas, siempre es necesario evaluar criticamente los resultados obtenidos.
