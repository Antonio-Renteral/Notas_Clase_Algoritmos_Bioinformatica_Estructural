En esta práctica trabajé con el archivo foldmason.pdb, que contiene varias estructuras superpuestas dentro de un solo archivo, y con foldmason_aa.fa,
que tiene el alineamiento de las secuencias. El objetivo era calcular el porcentaje de identidad y el RMSD entre algunas parejas de estructuras y resumir
los resultados en una tabla.

Usé el código original prog3.1.py, pero fue necesario hacer algunos ajustes para que funcionara correctamente con los datos. Primero entendí que no eran
varios archivos PDB separados, sino un solo archivo con múltiples modelos dentro, organizados con bloques MODEL y ENDMDL. El script ya estaba pensado para eso,
pero había que asegurarse de que leyera bien todos los modelos.

Después aparecieron errores como IndexError, que ocurrían porque el alineamiento tenía más posiciones que las coordenadas reales disponibles en la estructura.
No todos los residuos del alineamiento tienen coordenadas Cα en el PDB, así que el código intentaba acceder a posiciones que no existían. Para solucionarlo,
agregué validaciones para no salirme del rango de las listas, ignorar posiciones con gaps y verificar que realmente existieran coordenadas antes de usarlas
para el cálculo del RMSD.

También modifiqué el script para que generara automáticamente un archivo .csv con los resultados, incluyendo los nombres de los dominios, el porcentaje de identidad
y el RMSD. Con esas correcciones el programa pudo ejecutarse completo y producir la tabla final.

En general, esta práctica me ayudó a entender mejor la relación entre alineamiento de secuencias y comparación estructural, y a ver que en bioinformática
muchas veces los datos no coinciden perfectamente, por lo que es importante validar bien la información antes de hacer cálculos.
