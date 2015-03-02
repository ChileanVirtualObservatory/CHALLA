import math
import numpy as np
from scipy.optimize import leastsq
# Model: a,b,x0,y0,v0,dvx,dvy,Dx,Dy,Dv,phi

def _find_max(image):
   num=image.max()
   index=np.unravel_index(image.argmax(),image.shape)
   return (num,index)

def _clump2gauss((a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel))
   sphi2=math.pow(math.sin(phi),2)
   cphi2=math.pow(math.cos(phi),2)
   s2phi=math.sin(2*phi)
   Dalp2=math.pow(Dalp,2)
   Ddel2=math.pow(Ddel,2)
   Dv2=math.pow(Dv,2)
   La=cphi2/Dalp2 + sphi2/Ddel2 + math.pow(dvalp,2)/Dv2
   Ld=sphi2/Dalp2 + cphi2/Ddel2 + math.pow(dvdel,2)/Dv2
   Lb=-s2phi/(2*Dalp2) + s2phi/(2*Ddel2) + dvalp*dvdel/Dv2
   Lc=-dvalp/Dv2
   Le=-dvdel/Dv2
   Lf=1.0/Dv2
   L=np.array([[La,Lb,Lc],[Lb,Ld,Le],[Lc,Le,Lf]])
   mu=[alp0,del0,v0]
   return (a,b,mu,L)

def chi2(model, y, x):
   Gaussian(model,True)

def _modified_chi_leastsq(image,params,ymax,xmax):
   a=ymax
   b=0
   alp0=xmax[0]
   del0=xmax[1]
   v0=xmax[2]
   phi=0
   Dalp=params['beam_res']
   Ddel=params['beam_res']
   Dv=  params['velo_res']
   dvalp=0
   dvdel=0
   p0=(a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel)
   plsq = leastsq(chi2, p0, args=(y, x)) 
   print plsq[0]
   return plsq[0]


def gc_default_params():
   retval=dict()
   #TODO

def gc_gauss_clumps(ndd,params):
   C=[]
   stop=False
   image=ndd.data.copy();
   while not stop:
      (ymax,xmax)=_find_max(image)
      theta,image=_modified_chi_leastsq(image,params,ymax,xmax)
      C.append(theta)
      stop = (np.norm(image) < params['threshold'])
   return C
