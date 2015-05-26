import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *


def set_params():
   MaxJump=4
   MinDip=3*noise
   threshold=1


def create_caa(cube):
   caa=np.zeros_like(cube.data)
   # Check threshold
   return caa


def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)
       

def check_merge():
   return


def max_gradient(pixel):
    return 


def walkup(start,caa): 
   return 


def fellWallker(orig_cube):
   cube=cube.deepcopy(orig_cube)
   caa=create_caa(cube)
   data=cube.data
   shape=data.shape
   ()

   for i in range()

   
   return caa