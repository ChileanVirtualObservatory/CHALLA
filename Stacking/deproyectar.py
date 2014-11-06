#!/usr/bin/python

import math
import numpy as np
from info_imagen import distancia

def deproyectar(x,y,NAXIS1,NAXIS2,borde,matriz,radio):

    matriz_final = np.zeros(((y+(math.cos(90*math.pi/180))*radio)*2+2,radio*2+2))

    for i in range(NAXIS1):
        for j in range(NAXIS2):

            a = b = c = d = 0
            for r in range (len(borde)-1):

                if ( r % 2 != 0):
                    continue
                if i < borde[r]:
                    a = 1
                if i > borde[r]:
                    b = 10
                if j > borde[r+1]:
                    c = 100
                if j < borde[r+1]:
                    d = 1000
                if i == borde[r] and j == borde[r+1]:
                    a,b,c,d = 1,10,100,1000
                    break

            if (a + b + c + d != 1111 or x == i and y == j):
                continue

            angulo = math.asin(abs(j - abs(y))/distancia(x,y,i,j))*180/math.pi
            # Para mejor comprension, hay que leer como fila columna

            while(1):
                if i > x and j > y:
                    angulo = angulo
                    break
                if i < x and j > y:
                    angulo = 360 - angulo
                    break
                if i < x and j < y:
                    angulo = 180 + angulo
                    break
                if i > x and j < y:
                    angulo = 180 - angulo
                    break
                if i == x and j > y:
                    angulo = 0
                    break
                if i == x and j < y:
                    angulo = 180
                    break
                if i > x and j == y:
                    angulo = 90
                    break
                if i < x and j == y:
                    angulo = 270
                    break


            fil = x + (math.sin(angulo*math.pi/180))*radio
            col = y + (math.cos(angulo*math.pi/180))*radio

            matriz_final[int(fil),int(col)] = matriz[i][j]
            matriz_final[int(fil),math.ceil(col)] = matriz[i][j]
            matriz_final[math.ceil(fil),int(col)] = matriz[i][j]
            matriz_final[math.ceil(fil),math.ceil(col)] = matriz[i][j]

    return matriz_final