#!/usr/bin/python

import numpy as np

def alinear(x1,y1,x2,y2,dim_1_x,dim_1_y,dim_2_x,dim_2_y, actual): # 1: Referencia 2: Pos. de la imagen actual
    diferencia_x = x1 - x2
    diferencia_y = y1 - y2

    matriz_final = np.zeros((dim_1_x,dim_1_y))

    for i in range(dim_2_x):
        for j in range(dim_2_y):
            x = i + diferencia_x
            y = j + diferencia_y

            if ( x > 0 and x < dim_1_x and y > 0 and y < dim_1_y):
                matriz_final[x][y] = actual[i][j]

    return matriz_final