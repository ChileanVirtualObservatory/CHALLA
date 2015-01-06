#!/usr/bin/python

import math
from astropy.nddata import NDData
import numpy as np
from astropy.io import fits

def distancia(x1,y1,x2,y2):
    """
    mide la distancia entre 2 puntos utilizando teorema de pitagoras 

    Args:
        - x1: coordenada x del primer punto
        - y1: coordenada y del primer punto
        - x2: coordenada x del segundo punto
        - y2: coordenada y del segundo punto

    return: 
        - valor de la distancia calculada


    cambios: 
        - utilizar astropy.units para mejor representacion del dato
    """
    
    return math.sqrt(math.pow((x2-x1),2) + math.pow((y2-y1),2))

def centro_geometrico(dis,x1,y1,x2,y2,angulo):
    """
    calcula el centro geometrico de la elipse

    Args:
        - dis: longitud del eje mayor
        - x1,y1,x2,y2: coordenadas de los puntos extremos del eje mayor
        - angulo: angulo de inclinacion (DUDAAAAAA)

    return: 
        - punto_central = [x,y] 
        - punto de ubicacion del eje geometrico (x,y)
    """
    # Corresponde al centro geometrico de la figura (mitad eje mayor)

    if (angulo < 90):
        punto_central = [round(x1 + 0.5*math.cos(angulo*math.pi/180)*dis ), round(y1 + 0.5*math.sin(angulo*math.pi/180)*dis)]
    if (angulo > 90):
        punto_central = [round(x2 - 0.5*math.cos(angulo*math.pi/180)*dis ), round(y2 - 0.5*math.sin(angulo*math.pi/180)*dis)]
    if (angulo == 90):
        punto_central = [abs(x1-x2)/2,y1]

    return punto_central

def angulo_a_rotar(x1,y1,x2,y2):
    """
    calcula el angulo para que al rotar el eje mayor, este quede paralelo al eje x

    Args:
        - x1: coordenada x del primer punto del eje mayor
        - y1: coordenada y del primer punto del eje mayor
        - x2: coordenada x del segundo punto del eje mayor
        - y2: coordenada y del segundo punto del eje mayor

    return: 
        - angulo encontrado
    """
    # Indica el angulo a rotar para dejar el eje mayor paralelo al eje x
    if (((y2 - y1)*1.0) == 0):
        return 90
    else:
        return 180 - math.degrees(math.atan(abs(((y2 - y1)*1.0)/((x2 - x1)*1.0))))


def eje_mayor(borde):
    """
    encuentra los puntos separados a una mayor distancia en el borde de la elipse (eje mayor)

    Args:
        - borde: puntos correspondientes al borde de la elipse

    return: 
        - info_elipse = [Dimension,x1,y1,x2,y2] == arreglo con el valor numerico del eje mayor y los puntos extremos de este eje

    cambios: 
        - astropy.units
        - dato utilizado para representar "borde"


    """
    borde = NDData(borde) 
    info_elipse = np.zeros((4,))
    info_elipse = NDData(info_elipse)
    info_elipse.meta['distance'] = 0

    for i in range (0, len(borde.data)-2):
        if ( i % 2 != 0):
            continue
        for j in range (i+2, len(borde.data)-1):
            if ( j % 2 != 0):
                continue
            dis = distancia(borde.data[i],borde.data[i+1],borde.data[j],borde.data[j+1])
            if dis > info_elipse.meta['distance']:
                info_elipse.meta['distance'] = dis          # Distancia del eje mayor
                info_elipse.data[0] = borde.data[i]     # (X1,
                info_elipse.data[1] = borde.data[i+1]   #  Y1) Eje mayor
                info_elipse.data[2] = borde.data[j]     # (X2,
                info_elipse.data[3] = borde.data[j+1]   #  Y2) Eje mayor
    return info_elipse

if __name__ == "__main__":
    vec = [0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    info = eje_mayor(vec)



