#!/usr/bin/python

import numpy as np
import math
import info_imagen

def rotar(matriz, NAXIS1, NAXIS2, angulo):
    if (angulo > 360 or angulo < 1):
        print "<Error: Imagen no rotada, angulo no permitido>"
        return matriz
    # ------ PARA 0 NO ES NECESARIO ROTAR    ------ #
    if (angulo == 0 or angulo ==360):
        return matriz
    
    # ------ PARA 90, 180 y 270 ES UNA SIMPLE TRASLACION DE PUNTOS ------ #
    
    if (angulo == 90):
        matriz_final = np.zeros((NAXIS2,NAXIS1))
        for i in range(NAXIS1):
            for j in range(NAXIS2):
                matriz_final[NAXIS2 - j -1][i] = matriz[i][j]
        return matriz_final

    if (angulo == 180):
        matriz_final = np.zeros((NAXIS1,NAXIS2))
        for i in range(NAXIS1):
            for j in range(NAXIS2):
                matriz_final[NAXIS1 - i - 1][NAXIS2 - j -1] = matriz[i][j]
        return matriz_final

    if (angulo == 270):
        matriz_final = np.zeros((NAXIS2,NAXIS1))
        for i in range(NAXIS1):
            for j in range(NAXIS2):
                matriz_final[j][i] = matriz[i][j]
        return matriz_final

    else:
        
        coseno = math.cos((angulo*math.pi)/180)
        seno = math.sin((angulo*math.pi)/180)
        
        punto_central_x = int(round(NAXIS1/2))
        punto_central_y = int(round(NAXIS2/2))
        
    
        # --- Para rotar sobre el centro de la imagen, hay que hacer una pequena traslacion --- #
        # --- Conociendo la distancia del origen al centro de la imagen es suficiente       --- #
        distancia_centro = int(round(info_imagen.distancia(0,0,punto_central_x,punto_central_y))) - 1
        
        # --- PUNTO MAS NEGATIVO EN X Y EN Y ---------------------- #
        # --- ESTO ES PARA DEJAR TODAS LAS POSICIONES POSITIVAS --- #
        vec = [0,0,NAXIS1,NAXIS2,NAXIS1,0,0,NAXIS2]
        fila_mas_negativa = columna_mas_negativa = 0
        fila_mas_positiva = columna_mas_positiva = 0

        for i in range(7):
            alfa = (vec[i]-distancia_centro)*coseno - (vec[i+1]-distancia_centro)*seno
            beta = (vec[i]-distancia_centro)*seno + (vec[i+1]-distancia_centro)*coseno

            if (alfa < fila_mas_negativa):
                fila_mas_negativa = int(math.ceil(alfa))
            if (alfa > fila_mas_positiva):
                fila_mas_positiva = int(math.ceil(alfa))
            if (beta < columna_mas_negativa):
                columna_mas_negativa = int(math.ceil(beta))
            if (beta > columna_mas_positiva):
                columna_mas_positiva = int(math.ceil(beta))

        distancia_1 = fila_mas_positiva + abs(fila_mas_negativa)
        distancia_2 = columna_mas_positiva + abs(columna_mas_negativa)
        matriz_final = np.zeros((distancia_1+1,distancia_2+1))

        for x in range(NAXIS1):
            for y in range(NAXIS2):
            
                # ---- a X e Y hay que restarle y luego sumarle la traslacion -- #
                a = ((x-distancia_centro)*coseno - (y-distancia_centro)*seno ) + abs(fila_mas_negativa)
                b = ((x-distancia_centro)*seno + (y-distancia_centro)*coseno ) + abs(columna_mas_negativa)

                bandera_decimal_a = 100
                bandera_decimal_b = 100

                if( a - int(a) != 0):
                    bandera_decimal_a = 101
                if( b - int(b) != 0):
                    bandera_decimal_b = 110
            
                #Ya que en python no existe switch, se hace artesanalmente
                suma_banderas = bandera_decimal_a + bandera_decimal_b
                
                while(1):
                    
                    porcentaje_columna_derecha = porcentaje_columna_izquierda = 0
                    porcentaje_fila_abajo = porcentaje_fila_arriba = 0
                    
                    porcentaje_fila_arriba = abs(abs(a) - int(abs(a)))
                    porcentaje_fila_abajo  = 1 - porcentaje_fila_arriba
                    porcentaje_columna_derecha = abs(abs(b) - int(abs(b)))
                    porcentaje_columna_izquierda = 1 - porcentaje_columna_derecha

                    
                    #Solo A es decimal
                    if(suma_banderas == 201):
                        matriz_final[int(a)][b] += porcentaje_fila_abajo*matriz[x][y]
                        matriz_final[math.ceil(a)][b] += porcentaje_fila_arriba*matriz[x][y]
                        break
                                                                                    
                    #Solo B es decimal
                    if(suma_banderas == 210):
                        matriz_final[a][int(b)] += porcentaje_columna_izquierda*matriz[x][y]
                        matriz_final[a][math.ceil(b)] += porcentaje_columna_derecha*matriz[x][y]
                        break

                    #Ambos son decimales
                    if(suma_banderas == 211):
                        matriz_final[int(a)][int(b)] += porcentaje_fila_abajo*porcentaje_columna_izquierda*matriz[x][y]
                        matriz_final[math.ceil(a)][math.ceil(b)] += porcentaje_fila_arriba*porcentaje_columna_derecha*matriz[x][y]
                        matriz_final[int(a)][math.ceil(b)] += porcentaje_fila_abajo*porcentaje_columna_derecha*matriz[x][y]
                        matriz_final[math.ceil(a)][int(b)] +=  porcentaje_fila_arriba*porcentaje_columna_izquierda*matriz[x][y]
                        break
                    
                    #Ambos son enteros
                    if(suma_banderas == 200):
                        matriz_final[a][b] = matriz[x][y]
                        break
                            
        return matriz_final