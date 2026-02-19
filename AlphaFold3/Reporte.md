Durante este ejercicio trabaje en el modelado de un complejo proteina-ADN utilizando AlphaFold 3. En este caso seleccione el
factor de transcripcion LuxR, el cual funciona como homodimero y participa en regulacion genica dependiente de quorum sensing
en bacterias.

Elegi LuxR porque forma dimeros que reconocen secuencias especificas en el ADN conocidas como lux box. Esto lo hace adecuado para
modelar un complejo proteina-ADN con simetria parcial en el sitio de union.

Para el modelo utilice dos cadenas identicas de LuxR (cadena A y cadena B) para representar el homodimero. Como secuencia de ADN
emplee una lux box de 20 pares de bases:

5'- ACCTGTAGGATCGTACAGGT -3'  
3'- TGGACATCCTAGCATGTCCA -5'

Esta secuencia presenta una simetria parcial que permite que cada monomero del dimero interactue con una mitad del sitio de union.
No se incluyeron iones metalicos, ya que LuxR no requiere cofactores metalicos estructurales para unirse al ADN. Sin embargo,
es importante considerar que en condiciones biologicas reales LuxR se activa mediante la union a una molecula senal (AHL),
la cual estabiliza su estructura.

Una vez generada la prediccion en AlphaFold 3, evalue la calidad del modelo utilizando las metricas proporcionadas por el servidor.
Los valores de pLDDT fueron altos en el dominio de union al ADN, lo que sugiere una buena confianza estructural en la regio
funcional del complejo. Algunas regiones terminales mostraron menor confianza, lo cual es consistente con posibles segmentos
flexibles o parcialmente desordenados.

El mapa PAE mostro baja incertidumbre en la interfaz entre las dos cadenas de LuxR, indicando que el modelo predice correctamente
la dimerizacion. Tambien se observo baja incertidumbre en la region de contacto entre la proteina y el ADN, lo que sugiere que el
posicionamiento relativo es coherente.

Posteriormente el modelo fue analizado en SwissModel Assess para evaluar su calidad global. Los resultados indicaron una estructura
razonable, sin desviaciones geometricas importantes, lo que respalda la consistencia del complejo predicho.

En general, el ejercicio permitio observar como AlphaFold 3 puede modelar no solo estructuras proteicas individuales, sino tambien
complejos proteina-ADN. Ademas, fue interesante analizar como la simetria del sitio lux box favorece la union del homodimero y como
las metricas de confianza ayudan a interpretar la calidad estructural del modelo obtenido.
