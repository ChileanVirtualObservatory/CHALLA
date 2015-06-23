import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *


def get_params():
   maxJump=4
   minDip=3
   minSize=5
   flatSlope=0.1
   seaLevel=0.
   return {'maxJump':maxJump,'minDip':minDip,'minSize':minSize,'flatSlope':flatSlope,'seaLevel':seaLevel}


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
   """
   We will now examine the 3x3x3 cube of neighbouring pixels to determine
   where we step to next. We will choose the neighbour with the greatest
   gradient.
   """
   max_pos=pos
   max_grad=0.
   dist=[np.sqrt(1),np.sqrt(2),np.sqrt(3)]
   shape=data.shape

   for i in range(-1,2):
      for j in range(-1,2):
         for k in range(-1,2):
            neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point
            out1=neigh[0]<0 or neigh[1]<0 or neigh[2]<0 #out condition 1
            out2=neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2] #out condition 2
            d=dist[abs(i)+abs(j)+abs(k)-1] #distance to occupy

            if i==j==k==0:
               #dont check it again
               continue
            elif out1 or out2:
               #don't go outside of cube
               continue
            elif caa[neigh]==-1:
               #don't check unusable pixels
               continue
            else:
               grad=(data[neigh]-data[pos])/d #calculate gradient
               if grad>max_grad:
                  max_pos=neigh
                  max_grad=grad
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
            out1=neigh[0]<0 or neigh[1]<0 or neigh[2]<0 #out condition 1
            out2=neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2] #out condition 2

            if abs(i)<=1 and abs(j)<=1 and abs(k)<=1:
               #don't check it again
               continue
            elif out1 or out2:
               #don't go outside of cube
               continue
            elif caa[neigh]==-1:
               #don't check unusable pixels
               continue
            elif data[neigh]>max_val:
               max_pos=neigh
               max_val=data[neigh]
   return max_pos

"""
return dictionaries peaks and cols of the form:
peaks={clumpId:peak,...}
cols={clumpId:{neighId:col,...},...}
"""
def clump_structs(clump,data,caa):
   #Initialize peaks ans cols dictionaries
   peaks=dict()
   cols=dict()
   for clumpId in clump:
      cols[clumpId]=dict()
   shape=data.shape

   for clumpId,pixels in clump.items():
      peaks[clumpId]=0. #max value in current clump
      for pos in pixels:
         #First, find the peak of clump
         if data[pos]>peaks[clumpId]:
            peaks[clumpId]=data[pos]
         #Second, verify if it is a border clump pixel
         isBorder=False
         for i in range(-1,2):
            for j in range(-1,2):
               for k in range(-1,2):
                  neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point
                  out1=neigh[0]<0 or neigh[1]<0 or neigh[2]<0 #out condition 1
                  out2=neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2] #out condition 2
                  if i==j==k==0:
                     #don't check itself
                     continue
                  elif out1 or out2:
                     #don't go outside of cube
                     continue
                  elif caa[neigh]>0 and caa[neigh]!=clumpId:
                     #border clump pixel found!
                     isBorder=True
                     neighId=caa[neigh]
                     if neighId not in cols[clumpId]:
                        cols[clumpId][neighId]=data[pos]
                     elif data[pos]>cols[clumpId][neighId]:
                        cols[clumpId][neighId]=data[pos]
                     break
               if isBorder: break
            if isBorder: break
   return peaks,cols

def merge(clump,peaks,cols,caa):
   """
   Enter an iterative loop in which we join clumps that have a small dip
   between them. Quit when an interation fails to merge any more clumps.
   """
   merge=True
   while merge:
      #Don't iterate again unless some merge occur
      merge=False
      for clumpId in clump:
         """
         Check all the clumps that adjoin clump "clumpId", and find the one
         that is separated from "clumpId" by the smallest dip (i.e. has the
         highest "col" ).
         """
         topcol=0.
         neighId=-1
         for tmpId,col in cols[clumpId].items:
            if col>topcol:
               topcol=col
               neighId=tmpId

         """
         If the peak value in clumpId is less than "mindip" above the col
         between clumpId and neighId, then we merge the clumps. */
         """
         if neighId!=-1 and peaks[clumpId]<topcol+minDip:
            #Merge it!
            merge=True
            break #dictionaries can't change size while iterating it 
      if merge:
         """
         Tell the neighbours of clumpId that they are now neighbours of neighId
         instead. But only do this if they are not already neighbours of neighId.
         """
         #merging cols
         tmpCols=cols.pop(clumpId)
         tmpCols.update(cols[neighId])
         cols[neighId]=tmpCols
         for tmpId in cols:
            if clumpId in cols[tmpId]:
               cols[tmpId][neighId]=cols[tmpId].pop(clumpId)

         #delete clumpId from peaks and update new peak of merged clumps
         tmpPeak=peaks.pop(clumpId)
         peaks[neighId]=np.max(tmpPeak,peaks[neighId])

         #Update caa
         for pos in clump[clumpId]:
            caa[clumpId]=neighId

         #Merge pixels in clump dictionary
         clump[neighId]+=(clump.pop(clumpId))
         print "Merged clumps {0} and {1}".format(clumpId,neighId)

   return clump,peaks,cols,caa


def walkup(pos,path,pathv,data,caa):
   next_pos=max_gradient(pos,data,caa)

   if caa[next_pos]>=1:
      """ 
      If the next pixel on the route is already assigned to a
      clump, we now know what peak we are heading towards so use that clump
      index for the route so far, and abandon the rest of the walk.
      """
      path.append(next_pos)
      pathv.append(data[next_pos])
      return path,pathv

   elif next_pos!=pos:
      """
      Otherwise, move to the next pixel on the route.
      """
      path.append(next_pos)
      pathv.append(data[next_pos])
      return walkup(next_pos,path,pathv,data,caa)

   else:
      """
      If there is no upward route from the current pixel, (i.e. if it is a
      local maximum), we check the slightly more extended neighbourhood for
      a higher pixel. We do this because a local maximum could be just a
      noise spike. 
      """   
      local_max=next_pos
      new_max=verify_peak(local_max,data,caa)
      if local_max==new_max:
         """
         If the central pixel is the highest pixel in the more extended
         neighbourhood, we have reached a peak.
         """
         return path,pathv
      else:
         #Just a noise peak, keep walking up   
         return walkup(new_max,path,pathv,data,caa)

      """ 
      The walk has now finished because we have either reached a peak which
      is the highest point in its neighbourhood, or we have joined an earlier
      route and thus know which peak we would get to if we were to carry on.
      """

def verify_flat(path,pathv,caa,flatSlope,seaLevel):
   """ 
   Find the average gradient, if possible, over the last 4 steps, assuming
   each step is the same length. If it is above the user-supplied limit,
   indicate we have got the length of the initial flat section of the
   route. 
   """
   valid=-1 #valid flag, to indicate from what pixel path is valid


   #In case of small paths (not in cupid version)
   if 1<len(path)<4:
         if pathv[0]>=seaLevel:
            valid=0
         else:
            avg=(pathv[-1]-pathv[0])/(len(path)-1)
            if avg>=flatSlope:
               valid=0
   #In case of more than 4 steps in path
   else:
      for i in range(len(path)-3):
         if pathv[i]>=seaLevel:
            #valid from i-pixel
            valid=i
            break   
         else:
            avg=(pathv[i+3]-pathv[i])/3
            if avg>=flatSlope:
               #valid from i-pixel
               valid=i
               break

   #update path and flat 
   if valid==-1:
      flat=path
      flatv=pathv
      path=list()
      pathv=list()
   else:
      flat=path[0:valid]
      flatv=pathv[0:valid]
      path=path[valid::]
      pathv=pathv[valid::]

   return path,pathv,flat,flatv



def fellWalker(orig_cube):
   cube=copy.deepcopy(orig_cube)
   syn=copy.copy(orig_cube)
   syn.data=np.empty_like(cube.data)
   data=cube.data

   """
   Fill the supplied caa array with -1 for all pixels which are below the
   threshold, or are bad. Fill all other pixels with zero to indicate that
   the pixel is "usable but not yet checked".
   """
   caa=create_caa(data)

   """
   Use a cellular automata to remove small isolated groups of usable
   pixels from the above "ipa" array. This allocates new memory for the
   cleaned up array.
   """

   #Some constants
   rms=compute_rms(data)
   noise=1.*rms
   flatSlope=0.01

   seaLevel=2.*rms


   clump=dict()
   shape=data.shape
   top_id=0 #top clump id

   """
   Scan through the caa array, looking for usable pixels which have not
   yet been assigned to a clump (i.e. have a value of zero in caa).
   """
   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            # If unusable or already assigned, skip it
            if caa[i,j,k]!=0:
               continue
            path=list() # Ascent path pixels positions
            pathv=list() #Ascent path pixel values
            pos=(i,j,k)
            path.append(pos)
            pathv.append(data[pos])


            """
            We now start walking away from this pixel up hill (up the steepest
            gradient). We store the vector index of all pixels visited on this walk
            in the route array. Whilst walking, we initially keep a record of the 
            average gradient over the last 4 steps. We count how many steps we travel 
            before we encounter an average gradient above the gradient limit. 
            However, this initial flat section is only ignored if the walk starts from 
            "sea level". Walks which start at a higher altitude are used in their entirety.
            """
            path,pathv=walkup(pos,path,pathv,data,caa)
            path,pathv,flat,flatv=verify_flat(path,pathv,caa,flatSlope,seaLevel)

            """
            We now assign the clump index found above to all the pixels visited on
            he route, except for any low gradient section at the start of the route,
            which is set unusable (-1). We ignore walks that were entirely on the
            coastal plain (indicated by a value of -2 for clump_index).
            """
            
            #if not empty
            if flat:
               #set pixels as unusable
               for pos in flat:
                  caa[pos]=-1

            #if after that, path is empty then continue
            if not path:
               continue

            """
            If this peak has not already been assigned to a clump, increment the number 
            of clumps and assign it the new clump index.
            """
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

            #Print some info
            print "id={0}, path: {1}".format(top_id,path)

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
   """
   Amalgamate adjoining clumps if there is no significant dip between the
   clumps.
   """
   print "\n Merge Stage"

   """
   Get the minimum dip between two adjoining peaks necessary for the two
   peaks to be considered distinct.
   """
   minDip=get_params()['minDip']

   """
   Create a dictionary structure describing all the clumps boundaries in the
   supplied caa array.
   """
   clumpStruct=fwPixelStructs(clump,data,caa)


   """
   Enter an iterative loop in which we join clumps that have a small dip
   between them. Quit when an interation fails to merge any more clumps.
   """
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

   """
   Smooth the boundaries between the clumps. This cellular automata replaces
   each output pixels by the most commonly occuring value within a 3x3x3
   cube of input pixels centred on the output pixel. Repeat this process
   a number of times as given by configuration parameter CleanIter.
   """

   ####some statistics
   print "\n SOME USEFUL DATA"
   print "Number of clumps:",nclump
   for clumpId,pixels in clump.items():
      print "Clump {0} has {1} pixels".format(clumpId,len(pixels))

   return caa,clump