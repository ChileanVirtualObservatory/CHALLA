import numpy as np
#from statistics import Gaussian
from scipy.optimize import minimize
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
   return (a,b,mu,L)

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



def bubble_eval(features,(a,mu,sigmas)):
   C=np.empty_like(features
   C[0]=((features[0] - mu[0])**2)/(2.0*sigmas[0]**2)
   C[1]=((features[1] - mu[1])**2)/(2.0*sigmas[1]**2)
   C[2]=((features[2] - mu[2])**2)/(2.0*sigmas[2]**2)
   quad=C.sum(axis=0)
   retval=a*np.exp()
   return retval


def chi2(model,features,values,w,value_max,feature_max,params):
   sys.stdout.write('.')
   sys.stdout.flush()
   su=values - gauss_eval(features,to_gauss(model))
   nf=len(su) - 11;
   t1 = (np.square(su)*w).sum()/nf
   t2 = (np.exp(-su)).sum()/nf
   t3 = np.square((model[2]-feature_max[0])/params['beam_size']) + np.square((model[3]-feature_max[1])/params['beam_size']) + np.square((model[4]-feature_max[2])/params['spe_res'])
   t4 = np.square(model[0]+model[1] - value_max)
   return(t1 + params['s0']*t2 + params['sc']*t3 + params['sa']*t4)


def bubble(model,features,values,w,value_max,feature_max,params):
   return

   
def next_bubble(cube,syn,params):
   # Non-blocking plot 
   plt.ion() 
   plt.clf() 
   
   (value_max,feature_max) = cube.max()

   # Setting bubble parameters
   a=params['rms']
   x0=feature_max[0]
   y0=feature_max[1]
   v0=feature_max[2]
   mu=(x0,y0,v0)
   res_vect=np.array([params['beam_size'],params['beam_size'],params['spe_res']])
   sx=res_vect[0]
   sy=res_vect[1]
   sv=res_vect[2]
  
   # Compute the weight vector and the feature space
   w_sigmas=params['weight_deltas']*res_vect # several times the resolution
   (features,sc_index) = cube.feature_space(feature_max,2*w_sigmas)

   # Plot current cube
   plt.subplot(2, 3, 4)
   plt.imshow(cube.stack())
   rect=plt.Rectangle((sc_index[0],sc_index[2]),sc_index[1]-sc_index[0]+1,sc_index[3]-sc_index[2]+1,alpha=1, facecolor='none')
   plt.gca().add_patch(rect)
   
   # Plot current subcube
   plt.subplot(2, 3, 1)
   plt.imshow(cube.stack(sc_index))

   # Compute the bubble cube
   bubble_values=bubble_eval(features,(a,mu,(sx,sy,sv)))
   bubble_cube=cube_data_unravel(bubble_values,sc_index)
   
   # Remove bubble from the real cube 
   cube.add(-bubble_cube,sc_index)
   # Add bubble to the synthetic cube
   syn.add(bubble_cube,sc_index)
   
   # Plot current cube
   plt.subplot(2, 3, 5)
   plt.imshow(cube.stack())
   plt.gca().add_patch(rect)
   
   # Plot current subcube
   plt.subplot(2, 3, 2)
   plt.imshow(cube.stack(sc_index))
   
   # Plot bubble
   plt.subplot(2, 3, 3)
   plt.imshow(cube_data_stack(bubble_cube))

   # Plot synthetic
   plt.subplot(2, 3, 6)
   plt.imshow(syn.stack())
   plt.gca().add_patch(rect)

   plt.show() 
   plt.pause(0.01)

   # Return the max position
   return mu


def gauss_clumps_params():
   retval=dict()
   retval['threshold']=0.000001
   retval['few_deltas']=5
   retval['weight_deltas']=10
   retval['s0']=1.0
   retval['sc']=1.0
   retval['sa']=1.0
   retval['rms']=0.00002
   return retval


def get_bubbles(orig_cube,params):
   cube=copy.deepcopy(orig_cube)
   syn=copy.copy(orig_cube)
   syn.data=np.empty_like(cube.data)
   C=[]
   stop=False
   norm=cube.data.mean()
   print "Initial Sum", norm
   params['beam_size']=abs(float(cube.meta['BMIN']))
   params['spe_res']=abs(float(cube.meta['CDELT3']))
   params['rms']=compute_rms(cube.data)
   while not stop:
      mu=next_bubble(cube,syn,params)
      C.append(mu)
      norm=cube.data.mean()
      print "Sum", norm
      stop = (norm < params['threshold'])
   
   # SHOW RESULTS
   plt.clf() 
   plt.subplot(1, 3, 2)
   plt.imshow(cube.stack())
   plt.subplot(1, 3, 1)
   plt.imshow(orig_cube.stack())
   plt.subplot(1, 3, 3)
   plt.imshow(syn.stack())
   plt.show()
   plt.pause(100)
   print "RESULTS"
   return C


def get_gaussian():
   return 