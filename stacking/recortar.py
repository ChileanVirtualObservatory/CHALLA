#!/usr/bin/python

import numpy as np
from astropy.nddata import NDData

def recortar(matriz, borde):
    """
    extrae de la imagen el objeto a estudiar

    Args: 
        - matriz: imagen completa
        - borde: frontera del objeto de estudio

    return: 
        - matriz solo con el objeto de estudio

    cambios: 
        - ingreso de la matriz por input de archivo fits (astropy.io)
        - Borde. 
    """
    matriz = NDData(matriz)
    borde = NDData(borde)
    fila_mayor = 0
    fila_menor = borde.data[0]
    columna_mayor = 0
    columna_menor = borde.data[1]

    for i in range(len(borde.data)-1):

        if i % 2 != 0:
            continue
        if fila_mayor < borde.data[i]:
            fila_mayor = borde.data[i]
        if fila_menor > borde.data[i]:
            fila_menor = borde.data[i]
        if columna_mayor < borde.data[i+1]:
            columna_mayor = borde.data[i+1]
        if columna_menor > borde.data[i+1]:
            columna_menor = borde.data[i+1]

    for i in range(len(borde.data)-1):

        if i % 2 != 0:
            continue
        borde.data[i] = borde.data[i] - fila_menor
        borde.data[i+1] = borde.data[i+1] - columna_menor


    matriz_final = np.zeros((fila_mayor - fila_menor + 1, columna_mayor - columna_menor + 1))
    matriz_final = NDData(matriz_final)
    for i in range(fila_menor, fila_mayor+1):
        for j in range(columna_menor, columna_mayor+1):
            matriz_final.data[i - fila_menor][j - columna_menor] = matriz.data[i][j]

    return matriz_final
