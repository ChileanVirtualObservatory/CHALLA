import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *

def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)

def get_params():
   maxJump=4
   minDip=3
   minSize=5
   flatSlope=0.
   return {'maxJump':maxJump,'minDip':minDip,'minSize':minSize,'flatSlope':flatSlope}

def max_gradient(pos,data,caa):
   max_pos=pos
   max_val=data[pos]
   shape=data.shape

   for i in range(-1,2):
      for j in range(-1,2):
         for k in range(-1,2):
            neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point

            if i==j==k==0:
               #dont check it again
               continue
            elif neigh[0]<0 or neigh[1]<0 or neigh[2]<0 or neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2]:
               #don't go outside of cube
               continue
            elif caa[neigh]==-1:
               #don't check unusable pixels
               continue
            elif data[neigh]>max_val:
               max_pos=neigh
               max_val=data[neigh]
   return max_pos

def verify_peak(pos,data,caa):
   max_pos=pos
   max_val=data[pos]
   shape=data.shape
   maxJump=get_params()['maxJump']

   for i in range(-maxJump,maxJump+1):
      for j in range(-maxJump,maxJump+1):
         for k in range(-maxJump,maxJump+1):
            neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point

            if abs(i)<=1 and abs(j)<=1 and abs(k)<=1:
               #don't check it again
               continue
            elif neigh[0]<0 or neigh[1]<0 or neigh[2]<0 or neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2]:
               #don't go outside of cube
               continue
            elif caa[neigh]==-1:
               #don't check unusable pixels
               continue
            elif data[neigh]>max_val:
               max_pos=neigh
               max_val=data[neigh]
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
         return path
      else:
         #Just a noise peak, keep walking up   
         return walkup(new_max,path,data,caa)


def fwMain(data):
   """ The distance from a central pixel centre to the centre of a neighbouring
   pixel, indexed by the number of axes (0, 1, 2 or 3) on which the two
   pixels have zero displacement  (sqrt(3), sqrt(2), sqrt(1), 0). """
   dist=np.array([np.sqrt[3],np.sqrt[2],np.sqrt[1],0])

   """ Initialise the highest index value in caa. """
   top_id=0

   """ See if the user wants to dump the clump assignment array after a
   specified peak is reached for the first time. The supplied value is
   the clump index. """


   """ See if the user wants to dump the clump assignment array after the
   walk from a specified pixel has been completed. The supplied value is
   the zero-based vector index of the start pixel. """

   """ Get the lowest data value to be considered.  """
   rms=compute_rms(data)
   noise=2.*rms

   """ Get the largest jump (in pixels) which can be made from a local maximum
   to a higher maximum in the neighbourhood of the first maximum. """
   maxJump=get_params()['maxJump']

   """ Get the difference in data value between adjacent pixels which marks
   the start of a walk. Any initial section of the walk which has an
   average gradient (measured over 4 steps) less than this value is not
   assigned to any clump. """
   flatSlope=get_params()['flatSlope']

   """ Walks which start at low level are ignored until they achieve the
   gradient specified by FLATSLOPE. Walks which start at high level are
   used in their entirety. Store the data value which marks the break
   between high and low level. """
   seaLevel=noise+2.*rms

   """ Fill the supplied "ipa" array with -1 for all pixels which are below the
   threshold, or are bad. Fill all other pixels with zero to indicate that
   the pixel is "usable but not yet checked". """
   caa=np.zeros_like(data).astype(int)
   shape=caa.shape
   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            if data[i,j,k]<noise:
               caa[i,j,k]=-1
   print "CAA initialized successfully!\n"

   """ Use a cellular automata to remove small isolated groups of usable
   pixels from the above "ipa" array. This allocates new memory for the
   cleaned up array. """

   """ Pre-allocate some memory for the route array which holds the 1D vector
   indices of the pixels traversed in the route from the initial pixel to
   the clump peak. """
   path=list()

   """ Store some useful data. """
   clump=dict()

   """ Scan through the caa array, looking for usable pixels which have not
   yet been assigned to a clump (i.e. have a value of zero in "work"). """
   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            # If unusable or already assigned, skip it
            if caa[i,j,k]!=0:
               continue
            pos=(i,j,k)
            path.append(pos)

   """ We now start walking away from this pixel up hill (up the steepest
   gradient). We store the vector index of all pixels visited on this walk
   in the route array. Put the current pixel in it and initialise its length
   to one, and then loop until we have reached the peak. Whilst walking, we
   initially keep a record of the average gradient over the last 4 steps.
   We count how many steps we travel before we encounter an average gradient
   above the gradient limit. However, this initial flat section is only
   ignored if the walk starts from "sea level". Walks which start at a
   higher altitude are used in their entirety. """
            path=walkup(pos,path,data,caa)
            print "new path:",path

            if caa[path[-1]]>0:
               #Ascent path reach an existing path
               path_id=caa[path[-1]]
               for pos in path:
                  if pos==path[-1]:
                     #Already asigned to clump
                     continue
                  caa[pos]=path_id
                  clump[path_id].append(pos)
            else:
               #A new ascent path
               top_id+=1
               path_id=top_id
               clump[path_id]=list()
               for pos in path:
                  caa[pos]=path_id
                  clump[path_id].append(pos)
            print "top_id",top_id






   return caa,top_id


def fellWalker(orig_cube):
   cube=copy.deepcopy(orig_cube)
   syn=copy.copy(orig_cube)
   syn.data=np.empty_like(cube.data)
   data=cube.data

   caa,top_id=fwMain(data)