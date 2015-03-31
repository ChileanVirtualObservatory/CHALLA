import numpy as np
from statistics import Gaussian
from scipy.optimize import root
import copy
import matplotlib.pyplot as plt
import sys
from sympy import symbols,cos,sin,exp
# Model: a,b,x0,y0,v0,dvx,dvy,sx,sy,sv,phi



#def to_gauss((a,b,alp0,del0,v0,phi,sx,sy,sv,dvalp,dvdel)):
#   sphi2=np.square(np.sin(phi))
#   cphi2=np.square(np.cos(phi))
#   s2phi=np.sin(2*phi)
#   sx2=np.square(sx)
#   sy2=np.square(sy)
#   sv2=np.square(sv)
#   La=cphi2/sx2 + sphi2/sy2 + np.square(dvalp)/sv2
#   Ld=sphi2/sx2 + cphi2/sy2 + np.square(dvdel)/sv2
#   Lb=-s2phi/(2*sx2) + s2phi/(2*sy2) + dvalp*dvdel/sv2
#   Lc=-dvalp/sv2
#   Le=-dvdel/sv2
#   Lf=1.0/sv2
#   L=np.array([[Lf,Le,Lc],[Le,Ld,Lb],[Lc,Lb,La]])
#   mu=[v0,del0,alp0]
#   return (a,b,mu,L)

#def _chiJacobian(Xv,Yv,w,params):
#   #s0 = 1
#   #sc = 1
#   #sa = 1
#   nf = Xv.size - 11
#   ymax = max(Yv)
#   # Jacobian of Y (parameters)
#   Jy = np.array([diff(Y,a), diff(Y,b), diff(Y,x0), diff(Y,y0), diff(Y,v0), diff(Y,phi),
#           diff(Y,sx), diff(Y,sy), diff(Y,sv), diff(Y,dvx), diff(Y,dvy)])
#   chiJ = [0,0,0,0,0,0,0,0,0,0,0]
#
#   for i in range(Xv.size):
#           xv,yv,vv = Xv[i]
#           Yfit = Y.subs(x,xv).subs(y,yv).subs(v,vv)
#           for j in range(11):
#                   # Evaluates derivative respect to j-variable
#                   aux = Jy[j].subs(x,xv).subs(y,yv).subs(v,vv)
#                   chiJ[j] += -2*w[i]*(Y[i]-Yfit)*aux + s0*np.exp(Yfit-Y[i])*aux
#
#   # Sum extra terms for a and b derivatives
#   (a, b) = symbols('a b')
#   chiJ[0] += 2*sa*(a+b-Ymax)
#   chiJ[1] += 2*sa*(a+b-Ymax)
#   chiJ[2] += 2*sc*(x0-x)
#   chiJ[3] += 2*sc*(a+b-Ymax)
#   chiJ[4] += 2*sc*(a+b-Ymax)
#   return chiJ


#def _gc_chi2(model, Xv, Yv ,w,params):
#   #sys.stdout.write('.')
#   #sys.stdout.flush()
#   #wd=params['weight_deltas']*resv
#   
#   G=Gaussian(to_gauss(model),True)
#   t1=np.square(y-yi_fit)*w
#   t2=np.exp(yi_fit-y)
#   t3=np.square(model[2]-xmax[2])/np.square(resv[2]) + np.square(model[3]-xmax[1])/np.square(resv[1]) + np.square(model[4]-xmax[0])/np.square(resv[0])
#   t4=np.square(model[0]+model[1] - ymax)
#   s0=params['s0']
#   sc=params['sc']
#   sa=params['sa']
#   val=t1 + s0*t2 + sc*t3 + sa*t4
#   return val


class GaussClumpModel:
   def __init__(self,params):
      # Parameters which define a Gaussian clump
      self.shape = symbols('a b x0 y0 v0 phi sx sy sv dvx dvy')
      (a, b, x0, y0, v0, phi, sx, sy, sv, dvx, dvy) = self.shape
      # Variables of intensity Gaussian function
      (self.x, self.y, self.v) = symbols('x y v')

      # Precision Matrix
      self.Lambda = np.array([
         [cos(phi)**2/sx**2 + sin(phi)**2/sy**2 + dvx**2/sv**2, -sin(2*phi)/(2*sx**2) + sin(2*phi)/(2*sy**2) + dvx*dvy/sv**2, -dvx/sv**2],
         [-sin(2*phi)/(2*sx**2) + sin(2*phi)/(2*sy**2) + dvx*dvy/sv**2, sin(phi)**2/sx**2 + cos(phi)**2/sy**2 + dvy**2/sv**2, -dvy/sv**2],
         [-dvx/sv**2, -dvy/sv**2, 1/sv**2]])

      # x-mu vector
      u = np.array([self.x-x0, self.y-y0, self.v-v0])

      # Gaussian function
      self.model  = a * exp(-0.5*np.dot(u,np.dot(self.Lambda,u))) + b
      self.params = params

   def clump(self,features):
       #TODO
       pass

   def chi2(self,features,values,w):
       #TODO
       pass
   
   def jac(self,features,values,w):
       #TODO
       pass
   
   def eval(self,expr,shape):
      #(a,b,x0,y0,v0,phi,sx,sy,sv,dvx,dvy)=shape
      for i in range(len(shape)):
         expr=expr.subs(self.shape[i],shape[i])
      return expr

def next_clump(cube,syn,factory,params):
   # Non-blocking plot 
   plt.ion() 
   
   (value_max,feature_max) = cube.max()
   a=value_max
   b=0
   
   # Initial guess: position and orientation
   x0=feature_max[2]
   y0=feature_max[1]
   v0=feature_max[0]
   phi=0
   
   # Initial guess: variances
   sigmas=params['few_deltas']*params['res_vect'] # few times the resolution
   sx=sigmas[2]
   sy=sigmas[1]
   sv=sigmas[0]
 
   # Initial guess: redshift
   dvalp=0
   dvdel=0
  
   # Compute the weight vector and the feature space
   w_sigmas=params['weight_deltas']*params['res_vect'] # several times the resolution
   (features,(v_lb,v_ub,y_lb,y_ub,x_lb,x_ub)) = cube.feature_space(feature_max,2*w_sigmas) 
   w_shape=(1,0,x0,y0,v0,0,w_sigmas[2],w_sigmas[1],w_sigmas[0],0,0)
   sym=factory.clump(features)
   w=factory.eval(sym,w_shape)

   # Compile first guess
   guess=[a,b,x0,y0,v0,phi,sx,sy,sv,dvalp,dvdel]

   # Unravel the values
   lss=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1]
   values=lss.ravel()
   
   # Construct the chi2 and chi_func symbolic expressions  
   chi2=factory.chi2(features,values,w)
   chi2_func=lambda shape: factory.eval(chi2,shape)

   jaco=factory.jac(features,values,w)
   jaco_func=lambda shape: factory.eval(jaco,shape)

   # OPTIMIZE
   res = root(chi2_func,guess,method='lm',jac=jaco_func)
   print "clump =", res.x

   # Clump in the cube
   clu=factory.eval(sym,res.x)
   M=clu.reshape((v_ub-v_lb+1,y_ub-y_lb+1,x_ub-x_lb+1))
   
   # Matrices for displaying results (NOT ALGORITHMIC CODE)
   ma=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=0)
   spe=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=(1,2))
   prof=M.sum(axis=0)
   vmin=ma.min()
   vmax=ma.max()
   plt.clf() 
   plt.subplot(2, 3, 1)
   plt.imshow(ma,vmin=vmin,vmax=vmax)
   plt.subplot(2, 3, 3)
   plt.imshow(prof)
   plt.subplot(2, 3, 2)

   # Remove clump from the real cube 
   cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1] -= M
   
   # Add clump to the synthetic cube
   syn[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1] += M

   # SHOW results (NOT ALGORITHMIC CODE)
   spe2=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=(1,2))
   plt.imshow(cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=0),vmin=vmin,vmax=vmax)
   plt.subplot(2, 3, 6)
   plt.imshow(syn.sum(axis=0))
   plt.subplot(2, 3, 5)
   plt.imshow(cube.data.sum(axis=0))
   plt.gca().add_patch(plt.Rectangle((x_lb,y_lb),x_ub-x_lb+1,y_ub-y_lb+1,alpha=1, facecolor='none'))
   plt.subplot(2, 3, 4)
   plt.plot(spe,"b")
   plt.plot(spe2,"r")
   plt.show() 
   plt.pause(0.01)

   # Return the clump parameters
   return res.x


def gauss_clumps_params():
   retval=dict()
   retval['threshold']=0.000001
   retval['few_deltas']=5
   retval['weight_deltas']=40
   retval['s0']=1.0
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
   params['res_vect']=np.array([abs(float(cube.meta['CDELT3'])),abs(float(cube.meta['BMIN'])),abs(float(cube.meta['BMIN']))])
   factory=GaussClumpModel(params)
   while not stop:
      theta=next_clump(cube,syn,factory,params)
      C.append(theta)
      norm=cube.data.mean()
      print "norm", norm
      stop = (norm < params['threshold'])
   
   # SHOW RESULTS
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
