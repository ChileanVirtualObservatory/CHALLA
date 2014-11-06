#!/usr/bin/python

import numpy as np

def recortar(matriz, borde):
    fila_mayor = 0
    fila_menor = borde[0]
    columna_mayor = 0
    columna_menor = borde[1]

    for i in range(len(borde)-1):

        if i % 2 != 0:
            continue
        if fila_mayor < borde[i]:
            fila_mayor = borde[i]
        if fila_menor > borde[i]:
            fila_menor = borde[i]
        if columna_mayor < borde[i+1]:
            columna_mayor = borde[i+1]
        if columna_menor > borde[i+1]:
            columna_menor = borde[i+1]

    for i in range(len(borde)-1):

        if i % 2 != 0:
            continue
        borde[i] = borde[i] - fila_menor
        borde[i+1] = borde[i+1] - columna_menor


    matriz_final = np.zeros((fila_mayor - fila_menor + 1, columna_mayor - columna_menor + 1))

    for i in range(fila_menor, fila_mayor+1):
        for j in range(columna_menor, columna_mayor+1):
            matriz_final[i - fila_menor][j - columna_menor] = matriz[i][j]

    return matriz_final
