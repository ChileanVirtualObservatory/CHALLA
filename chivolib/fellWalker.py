import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *


def set_params():
   MaxJump=4
   MinDip=3*noise
   threshold=1


def create_caa(data):
   caa=np.zeros_like(data)
   # Check invalid pixels (below threshold)
   return caa


def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)
       

def check_merge():
   return


def max_gradient(pos,data,caa):
   max_pos=pos
   max_val=data[pos]
   shape=data.shape

   for i in range(-1,2):
      for j in range(-1,2):
         for k in range(-1,2):
            if i==j==k=0:
               # dont check it again
               continue
            else if pos[0]+i<0 | pos[1]+j<0 | pos[2]+k<0 | pos[0]+i>=shape[0] | pos[1]+j>=shape[1] | pos[2]+k>=shape[2]:
               # don't go outside of cube
               continue
            else if caa[i,j,k]==-1:
               # don't check unusable pixels
               continue
            else if data[i,j,k]>max_val:
               max_pos=(i,j,k)
               max_val=data[i,j,k]
   return max_pos

def verify_peak(pos,data,caa):
   max_pos=pos
   max_val=data[pos]
   shape=data.shape

   for i in range(-MaxJump,MaxJump+1):
      for j in range(-MaxJump,MaxJump+1):
         for k in range(-MaxJump,MaxJump+1):

   return max_pos

def walkup(pos,path,data,caa,flag):
   next_pos=max_gradient(pos,data,caa)

   if caa[next_pos]>=1:
      # Another ascent path reached
      flag=caa[next_pos]
      return None

   else if next_pos!=pos:
      # Keep walking up
      path.append(next_pos)
      walkup(next_pos,path,data,caa)

   else:
      # Local peak reached
      path.append(next_pos)
      local_max=next_pos
      new_max=verify_peak(local_max,data,caa):
      if local_max==new_max:
         # Significant peak reached
         flag+=1
         return
      else:
         # Just a noise peak, keep walking up
         path.append(new_max)
         walking(new_max,path,data,caa,flag)


def fellWalker(orig_cube):
   cube=cube.deepcopy(orig_cube)
   data=cube.data
   caa=create_caa(data)
   shape=data.shape
   top_id=0 # top clump id
   path=list()

   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            # If unusable or already assigned, skip it
            if caa[i,j,k]!=0:
               continue
            path=list() # Ascent path pixels positions
            pos=(i,j,k)
            path.append(pos)
            flag=walkup(pos,path,data,caa)
   return caa