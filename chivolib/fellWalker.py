import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *


def set_params():
   MaxJump=4
   MinDip=3
   threshold=1
   return (MaxJump,MinDip,threshold)


def create_caa(data):
   caa=np.zeros_like(data)
   #Check invalid pixels (below threshold)
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
            if i==j==k==0:
               #dont check it again
               continue
            elif pos[0]+i<0 | pos[1]+j<0 | pos[2]+k<0 | pos[0]+i>=shape[0] | pos[1]+j>=shape[1] | pos[2]+k>=shape[2]:
               #don't go outside of cube
               continue
            elif caa[i,j,k]==-1:
               #don't check unusable pixels
               continue
            elif data[i,j,k]>max_val:
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
            if abs(i)<=1 and abs(j)<=1 and abs(k)<=1:
               #don't check it again
               continue
            elif pos[0]+i<0 | pos[1]+j<0 | pos[2]+k<0 | pos[0]+i>=shape[0] | pos[1]+j>=shape[1] | pos[2]+k>=shape[2]:
               #don't go outside of cube
               continue
            elif caa[i,j,k]==-1:
               #don't check unusable pixels
               continue
            elif data[i,j,k]>max_val:
               max_pos=(i,j,k)
               max_val=data[i,j,k]
   return max_pos

def walkup(pos,path,data,caa):
   next_pos=max_gradient(pos,data,caa)

   if caa[next_pos]>=1:
      #Another ascent path reached
      path.append(next_pos)
      return path

   elif next_pos!=pos:
      #Keep walking up
      path.append(next_pos)
      return walkup(next_pos,path,data,caa)

   else:
      #Local peak reached
      local_max=next_pos
      new_max=verify_peak(local_max,data,caa)
      if local_max==new_max:
         #Significant peak reached
         path.append(local_max)
         return path
      else:
         #Just a noise peak, keep walking up   
         return walkup(new_max,path,data,caa)


def fellWalker(orig_cube):
   cube=cube.deepcopy(orig_cube)
   data=cube.data
   caa=create_caa(data)
   shape=data.shape
   top_id=0 #top clump id

   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            # If unusable or already assigned, skip it
            if caa[i,j,k]!=0:
               continue
            path=list() # Ascent path pixels positions
            pos=(i,j,k)
            path.append(pos)
            path=walkup(pos,path,data,caa)

            if caa[path[-1]]>0:
               #Ascent path reach an existing path
               path_id=caa[path[-1]]
               for pos in path:
                  caa[pos]=path_id
            else:
               #A new ascent path
               top_id+=1
               path_id=top_id
               for pos in path:
                  caa[pos]=path_id
   return caa