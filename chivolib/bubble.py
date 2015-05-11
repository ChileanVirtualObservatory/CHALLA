import numpy as np
from numpy.linalg import *
#from statistics import Gaussian
from scipy.cluster.vq import *
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *
import matplotlib.cm as cm
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

def gauss_eval(features,(a,b,mu,L)):
   C=np.empty_like(features)
   C[0]=features[0] - mu[0]
   C[1]=features[1] - mu[1]
   C[2]=features[2] - mu[2]
   V=C*(L.dot(C))
   quad=V.sum(axis=0)
   v=np.exp(-quad/2.0)

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
   retval=b + a*v*np.sqrt(det(L/(2*np.pi)))
   return retval

def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)

def bubble_clump(orig_cube):
   scale=3.0
   a_mult=10.0  
   thsh=2.0
   wsigmas=2.0
   cube=copy.deepcopy(orig_cube)
   syn=copy.copy(orig_cube)
   syn.data=np.empty_like(cube.data)
   stop=False
   #norm=cube.data[cube.data>0].mean()
   #print "Initial Mean", norm
   #bsize=abs(float(cube.meta['BMAJ']))
   ang_res=abs(float(cube.meta['CDELT1']))
   spe_res=abs(float(cube.meta['CDELT3']))
   #ssize=(bsize/ang_res)*spe_res
   rms=compute_rms(cube.data)
   print rms
   a=rms*a_mult
   th=rms*thsh
   # Window
   #sx=bsize*scale
   #sy=bsize*scale
   #sv=ssize*scale
   sx=ang_res*scale
   sy=ang_res*scale
   sv=spe_res*scale
   window=2.0*wsigmas*np.array([sx,sy,sv])
   # Compile Bubble
   res=np.array([a,0,0,0,0,0,sx,sy,sv,0,0])
   bubble=to_gauss(res)
   vect=np.empty((0,3))
   while not stop:
      (value_max,feature_max) = cube.max()
      bubble[2]=feature_max
      (features,sc_index) = cube.feature_space(feature_max,window) 
      val_fit=gauss_eval(features,bubble)
      fit_cube=cube_data_unravel(val_fit,sc_index)
      cube.add(-fit_cube,sc_index)
      syn.add(fit_cube,sc_index)
      vect=np.vstack((vect,feature_max))
      #norm=cube.data[cube.data>0].mean()
      #thsh=cube.data[cube.data<0].mean()
      print "vmax=", value_max , "(th=",th,")"
      stop = (value_max < th)
   
   # SHOW RESULTS
   L=bubble[3]
   Sig=inv(L)
   kmax=30
   error=np.empty(kmax)
   v0=np.abs(vect[:,0]).max()
   v1=np.abs(vect[:,1]).max()
   v2=np.abs(vect[:,2]).max()
   vv=np.empty_like(vect)
   vv[:,0]=vect[:,0]/v0
   vv[:,1]=vect[:,1]/v1
   vv[:,2]=vect[:,2]/v2
   for k in range(kmax):
      print "k",k+1
      (codebook,dist)=kmeans(vv,k+1)
      (clus,ddis)=vq(vv,codebook)
      e1=0
      e2=0
      for i in range(k+1):
         v=vect[clus==i]
         if (v.size==0):
            continue
         for mv in v:
             e1=e1+gauss_eval(v.T,(1.0,0,mv,L/2.0)).sum()
         ccv=np.cov(v.T)
         LL=inv(Sig+ccv)
         LL2=2*ccv
         dd=np.sqrt(det(2*np.pi*LL2))
         if np.isnan(dd):
            print LL2
            print det(2*np.pi*LL2)
         tv=np.array([v0,v1,v2])*codebook[i]
         e2=e2+np.square(gauss_eval(v.T,(1.0,0,tv,LL)).sum())*dd
      print e1,e2
      error[k]=np.square(a)*(e1 - e2)
   k=np.argmin(error)
   print "Final K=",k+1
   (codebook,dist)=kmeans(vv,k+1)
   (clus,ddis)=vq(vv,codebook)
   plt.clf() 
   zext=[orig_cube.ra_axis[0],orig_cube.ra_axis[-1],orig_cube.dec_axis[0],orig_cube.dec_axis[-1]]
   yext=[orig_cube.ra_axis[0],orig_cube.ra_axis[-1],orig_cube.nu_axis[0],orig_cube.nu_axis[-1]]
   xext=[orig_cube.dec_axis[0],orig_cube.dec_axis[-1],orig_cube.nu_axis[0],orig_cube.nu_axis[-1]]
   plt.subplot(4, 4, 1)
   plt.imshow(orig_cube.stack(),aspect='auto',origin='lower',extent=zext)
   plt.subplot(4, 4, 5)
   plt.imshow(orig_cube.stack(axis=1),aspect='auto',origin='lower',extent=yext)
   plt.subplot(4, 4, 9)
   plt.imshow(orig_cube.stack(axis=2),aspect='auto',origin='lower',extent=xext)
   plt.subplot(4, 4, 2)
   plt.imshow(cube.stack(),aspect='auto',origin='lower',extent=zext)
   plt.subplot(4, 4, 6)
   plt.imshow(cube.stack(axis=1),aspect='auto',origin='lower',extent=yext)
   plt.subplot(4, 4, 10)
   plt.imshow(cube.stack(axis=2),aspect='auto',origin='lower',extent=xext)
   plt.subplot(4, 4, 3)
   plt.imshow(syn.stack(),aspect='auto',origin='lower',extent=zext)
   plt.subplot(4, 4, 7)
   plt.imshow(syn.stack(axis=1),aspect='auto',origin='lower',extent=yext)
   plt.subplot(4, 4, 11)
   plt.imshow(syn.stack(axis=2),aspect='auto',origin='lower',extent=xext)
   plt.subplot(4, 4, 13)
   plt.plot(error)
   plt.subplot(4, 4, 4)
   plt.xlim(zext[0],zext[1])
   plt.ylim(zext[2],zext[3])
   colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
   for i in range(k+1):
      v=vect[clus==i]
      plt.scatter(v[:,0],v[:,1], color=next(colors))
   plt.subplot(4, 4, 8)
   plt.xlim(yext[0],yext[1])
   plt.ylim(yext[2],yext[3])
   colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
   for i in range(k+1):
      v=vect[clus==i]
      plt.scatter(v[:,0],v[:,2], color=next(colors))
   plt.subplot(4, 4, 12)
   plt.xlim(xext[0],xext[1])
   plt.ylim(xext[2],xext[3])
   colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
   for i in range(k+1):
      v=vect[clus==i]
      plt.scatter(v[:,1],v[:,2], color=next(colors))
   plt.show()
   plt.pause(100)
   return C
