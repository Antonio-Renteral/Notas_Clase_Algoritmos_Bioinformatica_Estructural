#!/usr/bin/env python

from __future__ import print_function
from math import sqrt
import SVD
import csv

""" prog3.1 modificado
Calcula % identidad y RMSD para parejas de estructuras
superpuestas obtenidas con FoldMason.
Lee:
 - foldmason.pdb (varios MODEL)
 - foldmason_aa.fa (alineamiento multiple)
Genera:
 - resultados_identidad_RMSD.csv
"""

__author__  = 'Modificado para practica CATH'

# ------------------------------------------------------------
# 1) FUNCIONES
# ------------------------------------------------------------

def lee_fasta(filename):
    secuencias = []
    with open(filename) as f:
        nombre = ""
        sec = ""
        for linea in f:
            linea = linea.strip()
            if linea.startswith(">"):
                if sec != "":
                    secuencias.append((nombre, sec))
                    sec = ""
                nombre = linea[1:]
            else:
                sec += linea
        if sec != "":
            secuencias.append((nombre, sec))
    return secuencias


def lee_modelos_PDB(filename):
    modelos = []
    modelo_actual = []
    res = ''
    prev_resID = ''

    with open(filename,'r') as pdbfile:
        for line in pdbfile:

            if line.startswith("MODEL"):
                modelo_actual = []
                res = ''
                prev_resID = ''
                continue

            if line.startswith("ENDMDL"):
                if res != '':
                    modelo_actual.append(res)
                modelos.append(modelo_actual)
                continue

            if line.startswith("ATOM"):

                resID = line[17:26]

                if resID != prev_resID:
                    if res != '':
                        modelo_actual.append(res)
                    res = line
                else:
                    res += line

                prev_resID = resID

    return modelos


def coords_alineadas(align1,coords1,align2,coords2):

    total1,total2 = -1,-1
    align_coords1,align_coords2 = [],[]
    length = len(align1)

    for r in range(0, length):
        res1 = align1[r:r+1]
        res2 = align2[r:r+1]

        if(res1 != '-'): total1+=1
        if(res2 != '-'): total2+=1

        if(res1 == '-' or res2 == '-'): continue

        align_coords1.append( extrae_coords_atomo(coords1[total1],' CA ') )
        align_coords2.append( extrae_coords_atomo(coords2[total2],' CA ') )

    return (align_coords1,align_coords2)


def extrae_coords_atomo(res,atomo_seleccion):

    atom_coords = []
    for atomo in res.split("\n"):
        if(atomo[12:16] == atomo_seleccion):
            atom_coords = [ float(atomo[30:38]),
                            float(atomo[38:46]),
                            float(atomo[46:54]) ]
    return atom_coords


def calcula_superposicion_SVD(pdbh1,pdbh2):

    def calcula_centro(coords):
        centro = [0,0,0]
        for coord in coords:
            for dim in range(0,3):
                centro[dim] += coord[dim]
        for dim in range(0,3):
            centro[dim] /= len(coords)
        return centro

    def calcula_coordenadas_centradas(coords,centro):
        return [[c[0]-centro[0],
                 c[1]-centro[1],
                 c[2]-centro[2]] for c in coords]

    def calcula_coordenadas_rotadas(coords,rotacion):
        rcoords = [0,0,0]
        for i in range(0,3):
            tmp = 0.0
            for j in range(0,3):
                tmp += coords[j] * rotacion[i][j]
            rcoords[i] = tmp
        return rcoords

    coords1 = pdbh1['align_coords']
    coords2 = pdbh2['align_coords']

    centro1 = calcula_centro(coords1)
    centro2 = calcula_centro(coords2)

    ccoords1 = calcula_coordenadas_centradas(coords1,centro1)
    ccoords2 = calcula_coordenadas_centradas(coords2,centro2)

    matriz = [[0,0,0],[0,0,0],[0,0,0]]
    peso = 1.0/len(ccoords1)

    for i in range(0,3):
        for j in range(0,3):
            tmp = 0.0
            for k in range(0,len(ccoords1)):
                tmp += ccoords1[k][i] * ccoords2[k][j] * peso
            matriz[i][j]=tmp

    [U, Sigma, V] = SVD.svd( matriz )

    rotacion = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(0,3):
        for j in range(0,3):
            rotacion[i][j]= U[j][0]*V[i][0] + \
                            U[j][1]*V[i][1] + \
                            U[j][2]*V[i][2]

    rmsd = 0.0
    for n in range(0,len(coords1)):
        coords1_rot = calcula_coordenadas_rotadas(ccoords1[n],rotacion)
        for i in range(0,3):
            desv = ccoords2[n][i]-coords1_rot[i]
            rmsd += desv*desv

    rmsd /= len(coords1)
    return sqrt(rmsd)


def calcula_identidad(align1,align2):

    matches = 0
    length = 0

    for a,b in zip(align1,align2):
        if a != "-" and b != "-":
            length += 1
            if a == b:
                matches += 1

    return 100.0 * matches / length


# ------------------------------------------------------------
# 2) PROGRAMA PRINCIPAL
# ------------------------------------------------------------

secuencias = lee_fasta("foldmason/foldmason_aa.fa")
modelos = lee_modelos_PDB("foldmason.pdb")

resultados = []

print("Dominio1\tDominio2\t%Identidad\tRMSD")

n = len(modelos)

for i in range(0,n):
    for j in range(i+1,n):

        nombre1, align1 = secuencias[i]
        nombre2, align2 = secuencias[j]

        pdb1 = {'coords':modelos[i]}
        pdb2 = {'coords':modelos[j]}

        (pdb1['align_coords'],pdb2['align_coords']) = \
            coords_alineadas(align1,modelos[i],
                             align2,modelos[j])

        rmsd = calcula_superposicion_SVD(pdb2,pdb1)
        identidad = calcula_identidad(align1,align2)

        print("%s\t%s\t%.2f\t%.2f" %
              (nombre1,nombre2,identidad,rmsd))

        resultados.append([nombre1,nombre2,
                           round(identidad,2),
                           round(rmsd,3)])

# ------------------------------------------------------------
# 3) GUARDAR RESULTADOS EN CSV
# ------------------------------------------------------------

with open("resultados_identidad_RMSD.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Dominio1","Dominio2",
                     "%Identidad","RMSD"])
    writer.writerows(resultados)

print("\nArchivo 'resultados_identidad_RMSD.csv' creado correctamente.")