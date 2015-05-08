import numpy as np
#from statistics import Gaussian
from scipy.optimize import fmin_bfgs,check_grad
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *
#from sympy import symbols,cos,sin,exp,lambdify,diff
# Model: a,b,x0,y0,v0,dvx,dvy,sx,sy,sv,phi


def to_gauss((a,b,x0,y0,v0,phi,sx,sy,sv,dvx,dvy)):
   sphi2=np.square(np.sin(phi))
   cphi2=np.square(np.cos(phi))
   s2phi=np.sin(2*phi)
   sx2=np.square(sx)
   sy2=np.square(sy)
   sv2=np.square(sv)
   La=cphi2/sx2 + sphi2/sy2 + np.square(dvx)/sv2
   Ld=sphi2/sx2 + cphi2/sy2 + np.square(dvy)/sv2
   Lb=-s2phi/(2*sx2) + s2phi/(2*sy2) + dvx*dvy/sv2
   Lc=-dvx/sv2
   Le=-dvy/sv2
   Lf=1.0/sv2
   L=np.array([[La,Lb,Lc],[Lb,Ld,Le],[Lc,Le,Lf]])
   mu=[x0,y0,v0]
   return [a,b,mu,L]

def gauss_eval(features,(a,b,mu,L),with_b=True):
   C=np.empty_like(features)
   C[0]=features[0] - mu[0]
   C[1]=features[1] - mu[1]
   C[2]=features[2] - mu[2]
   V=C*(L.dot(C))
   quad=V.sum(axis=0)
   
   # Standarize to 0-1 in an intelligent way to remove numerical precision problems
   #print quad
   quad=quad - quad.min()
   
   #plt.figure()
   #plt.plot(quad)
   #print(quad)
   v=np.exp(-quad/2.0)
   #plt.figure()
   #plt.plot(v)
   
   # 0-1 STandarization 
   #print v

   if (v.max()==0):
     print "ERROOOORRR quad"
     print quad
     print "v"
     print v
     print "Lambda"
     print L
     print "mu"
     print mu
     plt.figure()
     plt.plot(quad)
     print(quad)
     plt.figure()
     plt.plot(v)
   
   #plt.pause(100)
   #v=v/v.max()
   retval=a*v 
   if with_b:
      retval= retval + b
   return retval

def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)

def bubble_clump(orig_cube):
   count=1
   scale=3.0
   a_div=2.0
   wsigmas=2.0
   thsh=0.00005
   cube=copy.deepcopy(orig_cube)
   syn=copy.copy(orig_cube)
   syn.data=np.empty_like(cube.data)
   C=[]
   stop=False
   norm=cube.data.mean()
   print "Initial Mean", norm
   bsize=abs(float(cube.meta['BMAJ']))
   ang_res=abs(float(cube.meta['CDELT1']))
   spe_res=abs(float(cube.meta['CDELT3']))
   ssize=(bsize/ang_res)*spe_res
   rms=compute_rms(cube.data)
   print rms
   a=rms/a_div

   # Window
   res_vect=np.array([])
   sx=bsize*scale
   sy=bsize*scale
   sv=ssize*scale
   window=2.0*wsigmas*np.array([sx,sy,sv])
   # Compile Bubble
   res=np.array([a,0,0,0,0,0,sx,sy,sv,0,0])
   bubble=to_gauss(res)# Remove clump from the real cube 
   while not stop:
      #plt.ion() 
      #plt.clf() 
      (value_max,feature_max) = cube.max()
      bubble[2]=feature_max
      (features,sc_index) = cube.feature_space(feature_max,window) 
      val_fit=gauss_eval(features,bubble,False)
      fit_cube=cube_data_unravel(val_fit,sc_index)

      #Plot original cube
      #plt.subplot(2, 3, 4)
      #plt.imshow(orig_cube.stack())
      #rect=plt.Rectangle((sc_index[0],sc_index[2]),sc_index[1]-sc_index[0]+1,sc_index[3]-sc_index[2]+1,alpha=1, facecolor='none')
      #plt.gca().add_patch(rect)
   
      #Plot current subcube
      #plt.subplot(2, 3, 1)
      #plt.imshow(cube.stack(sc_index))

      cube.add(-fit_cube,sc_index)
      # Add clump to the synthetic cube
      syn.add(fit_cube,sc_index)
   
      #Plot current cube
      #plt.subplot(2, 3, 5)
      #plt.imshow(cube.stack())
      #plt.gca().add_patch(rect)
   
      #Plot current subcube
      #plt.subplot(2, 3, 2)
      #plt.imshow(cube.stack(sc_index))
   
      #Plot clump
      #plt.subplot(2, 3, 3)
      #plt.imshow(cube_data_stack(fit_cube))

      #Plot synthetic
      #plt.subplot(2, 3, 6)
      #plt.imshow(syn.stack())
      #plt.gca().add_patch(rect)
      #if count%10==0:
      #   plt.show() 
      #   plt.pause(0.00001)
      # Return the clump parameters
      
      theta=res.copy()
      theta[1:4]=feature_max
      C.append(theta)
      norm=cube.data.mean()
      print "Mean", norm
      stop = (norm < thsh)
      count+=1
   
   # SHOW RESULTS
   print len(C)
   print cube.data.shape
   plt.clf() 
   plt.subplot(1, 3, 2)
   plt.imshow(cube.stack())
   plt.subplot(1, 3, 1)
   plt.imshow(orig_cube.stack())
   plt.subplot(1, 3, 3)
   plt.imshow(syn.stack())
   plt.show()
   plt.pause(100)
   return C
