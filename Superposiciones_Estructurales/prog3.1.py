#!/usr/bin/env python

from __future__ import print_function
from math import sqrt
import SVD
import os
import csv

# ------------------------------------------------------------
# 1) LEE FASTA
# ------------------------------------------------------------

def lee_fasta(filename):

    secuencias = {}
    with open(filename) as f:
        nombre = ""
        sec = ""

        for linea in f:
            linea = linea.strip()

            if linea.startswith(">"):
                if sec != "":
                    secuencias[nombre] = sec
                    sec = ""
                nombre = linea[1:]
            else:
                sec += linea

        if sec != "":
            secuencias[nombre] = sec

    return secuencias


# ------------------------------------------------------------
# 2) LEE PDB INDIVIDUAL
# ------------------------------------------------------------

def lee_coordenadas_PDB(filename):

    coords = []
    with open(filename,'r') as pdbfile:

        res = ''
        prev_resID = ''

        for line in pdbfile:

            if line.startswith('TER'):
                break

            if not line.startswith('ATOM'):
                continue

            resID = line[17:26]

            if resID != prev_resID:
                if res != '':
                    coords.append(res)
                res = line
            else:
                res += line

            prev_resID = resID

        if res != '':
            coords.append(res)

    return coords


# ------------------------------------------------------------
# 3) FUNCIONES AUXILIARES
# ------------------------------------------------------------

def extrae_coords_atomo(res,atomo_seleccion):

    for atomo in res.split("\n"):
        if atomo[12:16] == atomo_seleccion:
            return [ float(atomo[30:38]),
                     float(atomo[38:46]),
                     float(atomo[46:54]) ]

    return []


def coords_alineadas(align1,coords1,align2,coords2):

    total1,total2 = -1,-1
    align_coords1,align_coords2 = [],[]
    length = len(align1)

    if length != len(align2):
        return ([],[])

    for r in range(length):

        res1 = align1[r]
        res2 = align2[r]

        if res1 != '-': total1 += 1
        if res2 != '-': total2 += 1

        if res1 == '-' or res2 == '-':
            continue

        if total1 >= len(coords1) or total2 >= len(coords2):
            continue

        c1 = extrae_coords_atomo(coords1[total1],' CA ')
        c2 = extrae_coords_atomo(coords2[total2],' CA ')

        if c1 != [] and c2 != []:
            align_coords1.append(c1)
            align_coords2.append(c2)

    return (align_coords1,align_coords2)


def calcula_identidad(align1,align2):

    matches = 0
    length = 0

    for a,b in zip(align1,align2):
        if a != "-" and b != "-":
            length += 1
            if a == b:
                matches += 1

    if length == 0:
        return 0.0

    return 100.0 * matches / length


def calcula_superposicion_SVD(pdbh1,pdbh2):

    coords1 = pdbh1['align_coords']
    coords2 = pdbh2['align_coords']

    if len(coords1) == 0:
        return None

    def calcula_centro(coords):
        centro = [0,0,0]
        for coord in coords:
            for dim in range(3):
                centro[dim] += coord[dim]
        for dim in range(3):
            centro[dim] /= len(coords)
        return centro

    def calcula_coordenadas_centradas(coords,centro):
        return [[c[0]-centro[0],
                 c[1]-centro[1],
                 c[2]-centro[2]] for c in coords]

    def calcula_coordenadas_rotadas(coords,rotacion):
        rcoords = [0,0,0]
        for i in range(3):
            tmp = 0.0
            for j in range(3):
                tmp += coords[j] * rotacion[i][j]
            rcoords[i] = tmp
        return rcoords

    centro1 = calcula_centro(coords1)
    centro2 = calcula_centro(coords2)

    ccoords1 = calcula_coordenadas_centradas(coords1,centro1)
    ccoords2 = calcula_coordenadas_centradas(coords2,centro2)

    matriz = [[0,0,0],[0,0,0],[0,0,0]]
    peso = 1.0/len(ccoords1)

    for i in range(3):
        for j in range(3):
            tmp = 0.0
            for k in range(len(ccoords1)):
                tmp += ccoords1[k][i] * ccoords2[k][j] * peso
            matriz[i][j] = tmp

    [U, Sigma, V] = SVD.svd(matriz)

    rotacion = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
        for j in range(3):
            rotacion[i][j] = U[j][0]*V[i][0] + \
                             U[j][1]*V[i][1] + \
                             U[j][2]*V[i][2]

    rmsd = 0.0
    for n in range(len(coords1)):
        coords1_rot = calcula_coordenadas_rotadas(ccoords1[n],rotacion)
        for i in range(3):
            desv = ccoords2[n][i]-coords1_rot[i]
            rmsd += desv*desv

    rmsd /= len(coords1)
    return sqrt(rmsd)


# ------------------------------------------------------------
# 4) PROGRAMA PRINCIPAL
# ------------------------------------------------------------

secuencias = lee_fasta("foldmason/foldmason_aa.fa")

pdbs = {}

# lee todos los .pdb que estén en la misma carpeta
for archivo in os.listdir("."):
    if archivo.endswith(".pdb"):
        nombre = archivo.replace(".pdb","")
        pdbs[nombre] = lee_coordenadas_PDB(archivo)

print("Total secuencias:", len(secuencias))
print("Total PDBs:", len(pdbs))

resultados = []
nombres = list(pdbs.keys())

print("Dominio1\tDominio2\t%Identidad\tRMSD")

for i in range(len(nombres)):
    for j in range(i+1, len(nombres)):

        nombre1 = nombres[i]
        nombre2 = nombres[j]

        if nombre1 not in secuencias or nombre2 not in secuencias:
            continue

        align1 = secuencias[nombre1]
        align2 = secuencias[nombre2]

        pdb1 = {'coords': pdbs[nombre1]}
        pdb2 = {'coords': pdbs[nombre2]}

        (pdb1['align_coords'], pdb2['align_coords']) = \
            coords_alineadas(align1, pdb1['coords'],
                             align2, pdb2['coords'])

        rmsd = calcula_superposicion_SVD(pdb2,pdb1)
        identidad = calcula_identidad(align1,align2)

        if rmsd is not None:
            print("%s\t%s\t%.2f\t%.2f" %
                  (nombre1,nombre2,identidad,rmsd))

            resultados.append([nombre1,nombre2,
                               round(identidad,2),
                               round(rmsd,3)])

with open("resultados_identidad_RMSD.csv","w",newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Dominio1","Dominio2",
                     "%Identidad","RMSD"])
    writer.writerows(resultados)

print("\nArchivo 'resultados_identidad_RMSD.csv' creado correctamente.")