import numpy as np
from statistics import Gaussian
from scipy.optimize import root
import copy
import matplotlib.pyplot as plt
import sys
# Model: a,b,x0,y0,v0,dvx,dvy,Dx,Dy,Dv,phi

def to_gauss((a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel)):
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

def _gc_chi2(model, y, X,resv,params,xmax,ymax):
   sys.stdout.write('.')
   sys.stdout.flush()
   wd=params['weight_deltas']*resv
   G=Gaussian(to_gauss(model),True)
   wmod=(1,0,model[2],model[3],model[4],0,wd[2],wd[1],wd[0],0,0)
   W=Gaussian(to_gauss(wmod),True)
   yi_fit=G.evaluate(X)
   wi=W.evaluate(X)
   t1=np.square(y-yi_fit)*wi
   t2=np.exp(yi_fit-y)
   t3=np.square(model[2]-xmax[2])/np.square(resv[2]) + np.square(model[3]-xmax[1])/np.square(resv[1]) + np.square(model[4]-xmax[0])/np.square(resv[0])
   t4=np.square(model[0]+model[1] - ymax)
   s0=params['s0']
   sc=params['sc']
   sa=params['sa']
   val=t1 + s0*t2 + sc*t3 + sa*t4
   return val

def _modified_chi_leastsq(cube,params,ymax,xmax,syn):
   plt.ion() 
   a=ymax
   b=0
   alp0=xmax[2]
   del0=xmax[1]
   v0=xmax[0]
   phi=0
   resv=np.array([abs(float(cube.meta['CDELT3'])),abs(float(cube.meta['BMIN'])),abs(float(cube.meta['BMIN']))])
   sigmas=params['few_deltas']*resv
   Dalp=sigmas[2]
   Ddel=sigmas[1]
   Dv=sigmas[0]
   dvalp=0
   dvdel=0
   (X,(n0,n1,d0,d1,r0,r1)) = cube.feature_space(xmax,2*params['weight_deltas']*resv) 
   p0=[a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel]
   print (n0,n1,d0,d1,r0,r1)
   print "p0 = ", p0
   lss=cube.data[n0:n1+1,d0:d1+1,r0:r1+1]
   #res= leastsq(_gc_chi2, p0, args=(lss.ravel(),X,resv,params,xmax,ymax)) 
   res = root(_gc_chi2,p0,method='lm',args=(lss.ravel(),X,resv,params,xmax,ymax))
   print "clump =", res[0]
   G=Gaussian(to_gauss(res.x),True)
   M=G.evaluate(X,False).reshape((n1-n0+1,d1-d0+1,r1-r0+1))
   ma=cube.data[n0:n1+1,d0:d1+1,r0:r1+1].sum(axis=0)
   spe=cube.data[n0:n1+1,d0:d1+1,r0:r1+1].sum(axis=(1,2))
   prof=M.sum(axis=0)
   vmin=ma.min()
   vmax=ma.max()
   plt.clf() 
   plt.subplot(2, 3, 1)
   plt.imshow(ma,vmin=vmin,vmax=vmax)
   plt.subplot(2, 3, 3)
   plt.imshow(prof)
   plt.subplot(2, 3, 2)
   cube.data[n0:n1+1,d0:d1+1,r0:r1+1] -= M
   syn[n0:n1+1,d0:d1+1,r0:r1+1] += M
   spe2=cube.data[n0:n1+1,d0:d1+1,r0:r1+1].sum(axis=(1,2))
   plt.imshow(cube.data[n0:n1+1,d0:d1+1,r0:r1+1].sum(axis=0),vmin=vmin,vmax=vmax)
   plt.subplot(2, 3, 6)
   plt.imshow(syn.sum(axis=0))
   plt.subplot(2, 3, 5)
   plt.imshow(cube.data.sum(axis=0))
   plt.gca().add_patch(plt.Rectangle((r0,d0),r1-r0+1,d1-d0+1,alpha=1, facecolor='none'))
   plt.subplot(2, 3, 4)
   plt.plot(spe,"b")
   plt.plot(spe2,"r")
   plt.show() 
   plt.pause(0.01)
   return res.x


def gc_default_params():
   retval=dict()
   retval['threshold']=0.000001
   retval['few_deltas']=5
   retval['weight_deltas']=40
   retval['s0']=100
   retval['sc']=1.0
   retval['sa']=1.0

   return retval

def gauss_clumps(orig_cube,params):
   cube=copy.deepcopy(orig_cube)
   syn=np.empty_like(cube.data)
   C=[]
   stop=False
   norm=cube.data.mean()
   print "Initial Norm", norm
   while not stop:
      (ymax,xmax)= cube.max()
      print "ymax,xmax = ",ymax,xmax
      theta=_modified_chi_leastsq(cube,params,ymax,xmax,syn)
      C.append(theta)
      norm=cube.data.mean()
      print "norm", norm
      stop = (norm < params['threshold'])
   
   plt.clf() 
   plt.subplot(1, 3, 2)
   plt.imshow(cube.data.sum(axis=0))
   plt.subplot(1, 3, 1)
   plt.imshow(orig_cube.data.sum(axis=0))
   plt.subplot(1, 3, 3)
   plt.imshow(syn.sum(axis=0))
   plt.show()
   plt.pause(100)
   print "RESULTS"
   print C
   return C
