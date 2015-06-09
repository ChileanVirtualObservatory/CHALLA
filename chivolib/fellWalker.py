import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *


def get_params():
   maxJump=4
   minDip=3
   minSize=5
   return {'maxJump':maxJump,'minDip':minDip,'minSize':minSize}


def create_caa(data):
   caa=np.zeros_like(data).astype(int)
   #Check invalid pixels (below threshold)
   shape=caa.shape
   rms=compute_rms(data)
   threshold=3*rms # here for now

   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            if data[i,j,k]<threshold:
               caa[i,j,k]=-1
   print "CAA initialized successfully"
   return caa


def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)


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

def check_merge(clumpId,clump,data,caa):
   shape=data.shape
   max_clump=0 #max value in clump
   max_border=0 #max value in border
   neighId=-1
   for pos in clump[clumpId]:
      #First, find the peak of clump
      if data[pos]>max_clump:
         max_clump=data[pos]
      #Second, verify if it is a border clump pixel
      for i in range(-1,2):
         for j in range(-1,2):
            for k in range(-1,2):
               neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point
               if i==j==k==0:
                  #don't check itself
                  continue
               elif neigh[0]<0 or neigh[1]<0 or neigh[2]<0 or neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2]:
                  #don't go outside of cube
                  continue
               elif caa[neigh]>0 and caa[neigh]!=clumpId:
                  #border clump pixel found!
                  #Verify if it is the highest
                  border_val=(data[pos]+data[neigh])/2. #value at border
                  if border_val>max_border:
                     max_border=border_val
                     neighId=caa[neigh]
   return (max_clump-max_border,neighId)


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


def fellWalker(orig_cube):
   cube=copy.deepcopy(orig_cube)
   data=cube.data
   caa=create_caa(data)
   clump=dict()
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

   ###refine 1, removing small clumps
   print "\n Removing small clumps stage"
   minSize=get_params()['minSize']
   deleted=list() #deleted id's

   for clumpId,pixels in clump.items():
      if len(pixels)<=minSize:
         #Set pixels as unusable
         for pos in pixels:
            caa[pos]=-1
         del clump[clumpId]
         deleted.append(clumpId)
         print "removed clump {0}".format(clumpId)
      elif deleted:
         #If deleted have any index
         clump[deleted[0]]=clump.pop(clumpId)
         #update caa
         for pos in clump[deleted[0]]:
            caa[pos]=deleted[0]
         del deleted[0]
         deleted.append(clumpId)

   
   ####refine 2, merging clumps
   print "\n Merge Stage"
   minDip=get_params()['minDip']
   merge=True
   while merge:
      #Don't iterate again unless some merge occur
      merge=False
      for clumpId in clump:
         dip,neighId=check_merge(clumpId,clump,data,caa)
         if neighId!=-1 and dip<=minDip:
            #Merge it!
            merge=True
            break #dictionaries can't change size while iterating it 
      if merge:
         #Update caa
         for pos in clump[neighId]:
            caa[pos]=clumpId
         #Merge pixels
         clump[clumpId]+=(clump.pop(neighId))
         print "Merged clumps {0} and {1}".format(clumpId,neighId)
   
   #ordering clumpId's (sequential id's)
   seqId=1
   for clumpId in clump:
      if clumpId!=seqId:
         clump[seqId]=clump.pop(clumpId)
         #update caa
         for pos in clump[seqId]:
            caa[pos]=seqId
      seqId+=1

   ####some statistics
   print "\n SOME USEFUL DATA"
   print "Number of clumps:",len(clump)
   for clumpId,pixels in clump.items():
      print "Clump {0} has {1} pixels".format(clumpId,len(pixels))

   return caa,clump