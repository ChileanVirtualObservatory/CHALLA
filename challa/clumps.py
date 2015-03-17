import numpy as np
from functions import Gaussian
from scipy.optimize import leastsq
import copy
import matplotlib.pyplot as plt

# Model: a,b,x0,y0,v0,dvx,dvy,Dx,Dy,Dv,phi

def clump2gauss((a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel)):
   sphi2=np.square(np.sin(phi))
   cphi2=np.square(np.cos(phi))
   s2phi=np.sin(2*phi)
   Dalp2=np.square(Dalp)
   Ddel2=np.square(Ddel)
   Dv2=np.square(Dv)
   La=cphi2/Dalp2 + sphi2/Ddel2 + np.square(dvalp)/Dv2
   Ld=sphi2/Dalp2 + cphi2/Ddel2 + np.square(dvdel)/Dv2
   Lb=-s2phi/(2*Dalp2) + s2phi/(2*Ddel2) + dvalp*dvdel/Dv2
   Lc=-dvalp/Dv2
   Le=-dvdel/Dv2
   Lf=1.0/Dv2
   L=np.array([[Lf,Le,Lc],[Le,Ld,Lb],[Lc,Lb,La]])
   mu=[v0,del0,alp0]
   return (a,b,mu,L)

def gc_chi2(model, y, X,resv,params,xmax,ymax):
   print model
   wd=params['weight_deltas']*resv
   G=Gaussian(clump2gauss(model),True)
   wmod=(1,0,model[2],model[3],model[4],0,wd[2],wd[1],wd[0],0,0)
   W=Gaussian(clump2gauss(wmod),True)
   yi_fit=G.evaluate(X)
   #print len(yi_fit)
   #plt.plot(y[6100000:],"r")
   #plt.plot(yi_fit[6100000:],"b")
   #plt.show()
   #print "model"
   #print model
   #print "yifit"
   #print yi_fit
   wi=W.evaluate(X)
   #print y.shape, yi_fit.shape
   #print X.shape
   #print y.shape,yi_fit.shape
   t1=np.power(y-yi_fit,2)*wi
   t2=np.exp(yi_fit-y)
   t3=np.square(model[2]-xmax[2])/np.square(resv[2]) + np.square(model[3]-xmax[1])/np.square(resv[1]) + np.square(model[4]-xmax[0])/np.square(resv[0])
   t4=np.square(model[0]+model[1] - ymax)
   s0=params['s0']
   sc=params['sc']
   sa=params['sa']
   val=t1 + s0*t2 + sc*t3 + sa*t4
   return np.sqrt(val)
   #print np.sqrt(val.sum())
   #return rv

def _modified_chi_leastsq(cube,params,ymax,xmax):
   a=ymax
   b=0
   alp0=xmax[2]
   del0=xmax[1]
   v0=xmax[0]
   phi=0
   resv=np.array([10*float(cube.meta['CDELT3']),float(cube.meta['BMAJ']),float(cube.meta['BMAJ'])])
   sigmas=params['few_deltas']*resv
   Dalp=sigmas[2]
   Ddel=sigmas[1]
   Dv=sigmas[0]
   dvalp=0
   dvdel=0
   (X,(n0,n1,d0,d1,r0,r1)) = cube.feature_space(xmax,2*params['weight_deltas']*resv) 
   p0=[a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel]
   #print "leastsq params"
   #print p0
   #print y
   #print X
   print "LEAST!"
   #print cube.data.shape
   lss=cube.data[n0:n1+1,d0:d1+1,r0:r1+1]
   #print lss.shape
   res= leastsq(gc_chi2, p0, args=(lss.ravel(),X,resv,params,xmax,ymax)) 
   #print res[0]
   G=Gaussian(clump2gauss(res[0]),True)
   M=G.evaluate(X,False).reshape((n1-n0+1,d1-d0+1,r1-r0+1))
   #print (n0,n1,d0,d1,r0,r1)
   #print "M"
   #clum1=M.sum(axis=1)
   #plt.imshow(clum1)
   #plt.show()
   cube.data[n0:n1+1,d0:d1+1,r0:r1+1] -= M
   return res[0]


def gc_default_params():
   retval=dict()
   retval['threshold']=0.01
   retval['few_deltas']=30
   retval['weight_deltas']=100
   retval['s0']=0.3
   retval['sc']=0.3
   retval['sa']=0.3

   return retval

def gauss_clumps(orig_cube,params):
   cube=copy.deepcopy(orig_cube)
   C=[]
   stop=False
   while not stop:
      (ymax,xmax)= cube.max()
      print "ymax,xmax"
      print ymax,xmax
      theta=_modified_chi_leastsq(cube,params,ymax,xmax)

      C.append(theta)
      norm=np.linalg.norm(cube.data)
      print "norm"
      print norm
      stop = (norm < params['threshold'])
      plt.imshow(cube.data.sum(axis=0))
      plt.show()
   return C
