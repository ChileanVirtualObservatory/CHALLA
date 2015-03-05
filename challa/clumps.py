import math
import numpy as np
from scipy.optimize import leastsq
# Model: a,b,x0,y0,v0,dvx,dvy,Dx,Dy,Dv,phi

def clump2gauss((a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel))
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

def gc_chi2(model, y, X,resv,params,xmax,ymax):
   wd=params['weight_deltas']*resv
   G=Gaussian(clump2gauss(model),True)
   wmod=(1,0,model[2],model[3],model[4],0,wd[0],wd[1],wd[2],0,0)
   W=Gaussian(clump2gauss(wmod,True)
   yi_fit=G.evaluate(X)
   wi=W.evaluate(X)
   t1=np.power(y-yi_fit,2)*wi
   t2=np.exp(yi_fit-y)
   t3=math.pow(model[2]-xmax[0],2)/math.pow(resv[0],2) + math.pow(model[3]-xmax[1],2)/math.pow(resv[1],2) + math.pow(model[4]-xmax[2],2)/math.pow(resv[2],2)
   t4=math.pow(model[0]+model[1] - ymax,2)
   s0=params['s0']
   sc=params['sc']
   sa=params['sa']
   val=t1 + s0*t2 + sc*t3 + sa*t4
   print val
   return np.sqrt(val)

def _modified_chi_leastsq(cube,params,ymax,xmax):
   a=ymax
   b=0
   alp0=xmax[0]
   del0=xmax[1]
   v0=xmax[2]
   phi=0
   resv=[float(cube.meta['BMAJ']),float(cube.meta['BMAJ'],float(cube.meta['CDELT3']
   sigmas=params['few_deltas']*resv
   Dalp=sigmas[0]
   #params['few_deltas']*float(cube.meta['BMAJ'])
   Ddel=sigmas[1]
   #params['few_deltas']*float(cube.meta['BMAJ'])
   Dv=sigmas[2]
   #params['few_deltas']*float(cube.meta['CDELT3'])
   dvalp=0
   dvdel=0
   #wd=params['weight_deltas']*resv
   y,X = cube.subcube(xmax,2*params['weight_deltas']*resv)
   p0=(a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel)
   plsq = leastsq(gc_chi2, p0, args=(y, X,resv,params,xmax,ymax)) 
   print plsq[0]
   return plsq[0]


def gc_default_params():
   retval=dict()
   retval['threshold']=0.0000001
   retval['few_deltas']=3
   retval['weight_deltas']=50
   return retval

def gauss_clumps(orig_cube,params):
   cube=copy.deepcopy(orig_cube)
   C=[]
   stop=False
   while not stop:
      (ymax,xmax)= cube.max()
      (theta,image)=_modified_chi_leastsq(cube,params,ymax,xmax)
      C.append(theta)
      stop = (np.norm(image) < params['threshold'])
   return C
