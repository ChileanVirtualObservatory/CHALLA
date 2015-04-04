import numpy as np
#from statistics import Gaussian
from scipy.optimize import root
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
   v=np.exp(-quad/2)
   retval=a*v 
   if with_b:
      retval= retval + b
   return retval

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


#class GaussClumpModel:
#   def __init__(self,params):
#      # Parameters which define a Gaussian clump
#      self.shape = symbols('a b x0 y0 v0 phi sx sy sv dvx dvy')
#      (a, b, x0, y0, v0, phi, sx, sy, sv, dvx, dvy) = self.shape
#      # Variables of intensity Gaussian function
#      (self.x, self.y, self.v) = symbols('x y v')
#
#      # Precision Matrix
#      self.Lambda = np.array([
#         [cos(phi)**2/sx**2 + sin(phi)**2/sy**2 + dvx**2/sv**2, -sin(2*phi)/(2*sx**2) + sin(2*phi)/(2*sy**2) + dvx*dvy/sv**2, -dvx/sv**2],
#         [-sin(2*phi)/(2*sx**2) + sin(2*phi)/(2*sy**2) + dvx*dvy/sv**2, sin(phi)**2/sx**2 + cos(phi)**2/sy**2 + dvy**2/sv**2, -dvy/sv**2],
#         [-dvx/sv**2, -dvy/sv**2, 1/sv**2]])
#
#      # x-mu vector
#      u = np.array([self.x-x0, self.y-y0, self.v-v0])
#
#      # Gaussian function
#      self.model  = a * exp(-0.5*np.dot(u,np.dot(self.Lambda,u))) + b
#      M=self.model
#      self.params = params
#      self.fitfunc = np.vectorize(lambdify((self.x,self.y,self.v,a, b, x0, y0, v0, phi, sx, sy, sv, dvx, dvy),self.model))
#      self.jacomodel = np.array([diff(M,a), diff(M,b), diff(M,x0), diff(M,y0), diff(M,v0), diff(M,phi),
#           diff(M,sx), diff(M,sy), diff(M,sv), diff(M,dvx), diff(M,dvy)])
#      
#   def clump_func(self,features):
#       ff=lambda par: self.fitfunc(features[0],features[1],features[2],par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10])
#       return ff
#
def chi2(model,features,values,w,value_max,feature_max,params):
   sys.stdout.write('.')
   sys.stdout.flush()
   su=values - gauss_eval(features,to_gauss(model))
   t1 = np.square(su)*w
   t2 = np.exp(-su)
   t3 = np.square((model[2]-feature_max[0])/params['beam_size']) + np.square((model[3]-feature_max[1])/params['beam_size']) + np.square((model[4]-feature_max[2])/params['spe_res'])
   t4 = np.square(model[0]+model[1] - value_max)
   return(t1 + params['s0']*t2 + params['sc']*t3 + params['sa']*t4)

def next_clump(cube,syn,params):
   # Non-blocking plot 
   plt.ion() 
   plt.clf() 
   
   (value_max,feature_max) = cube.max()
   a=value_max
   b=0
   
   # Initial guess: position and orientation
   x0=feature_max[0]
   y0=feature_max[1]
   v0=feature_max[2]
   phi=0
   
   # Initial guess: variances
   res_vect=np.array([params['beam_size'],params['beam_size'],params['spe_res']])
   sigmas=params['few_deltas']*res_vect # few times the resolution
   sx=sigmas[1]
   sy=sigmas[1]
   sv=sigmas[2]
 
   # Initial guess: redshift
   dvalp=0
   dvdel=0
  
   # Compute the weight vector and the feature space
   w_sigmas=params['weight_deltas']*res_vect # several times the resolution
   (features,sc_index) = cube.feature_space(feature_max,2*w_sigmas) 
   w_shape=(1,0,x0,y0,v0,0,w_sigmas[2],w_sigmas[1],w_sigmas[0],0,0)
   w=gauss_eval(features,to_gauss(w_shape),False)
 
   #Plot current cube
   plt.subplot(2, 3, 4)
   plt.imshow(cube.stack())
   rect=plt.Rectangle((sc_index[0],sc_index[2]),sc_index[1]-sc_index[0]+1,sc_index[3]-sc_index[2]+1,alpha=1, facecolor='none')
   plt.gca().add_patch(rect)
   
   #Plot current subcube
   plt.subplot(2, 3, 1)
   plt.imshow(cube.stack(sc_index))

   # Compile first guess
   guess=[a,b,x0,y0,v0,phi,sx,sy,sv,dvalp,dvdel]

   # Ravel the values
   values=cube.ravel(sc_index)
   
   # Pack all args of chi2
   chi2_args=(features,values,w,value_max,feature_max,params)
   # Construct the chi2 and chi_func functions

   # OPTIMIZE
   #res = root(chi2_func,guess,method='lm',jac=jaco_func,args=chi2_args)
   res = root(chi2,guess,method='lm',args=chi2_args)
   print "clump =", res.x

   # Clump values
   val_fit=gauss_eval(features,to_gauss(res.x),False)
   fit_cube=cube_data_unravel(val_fit,sc_index)
   # Remove clump from the real cube 
   cube.add(-fit_cube,sc_index)
   # Add clump to the synthetic cube
   syn.add(fit_cube,sc_index)
   
   #Plot current cube
   plt.subplot(2, 3, 5)
   plt.imshow(cube.stack())
   plt.gca().add_patch(rect)
   
   #Plot current subcube
   plt.subplot(2, 3, 2)
   plt.imshow(cube.stack(sc_index))
   
   #Plot clump
   plt.subplot(2, 3, 3)
   plt.imshow(cube_data_stack(fit_cube))

   #Plot synthetic
   plt.subplot(2, 3, 6)
   plt.imshow(syn.stack())
   plt.gca().add_patch(rect)

   #END
   #M=clu.reshape((v_ub-v_lb+1,y_ub-y_lb+1,x_ub-x_lb+1))
   
   #cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1] -= M
   #syn[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1] += M
   
   # Matrices for displaying results (NOT ALGORITHMIC CODE)
   #ma=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=0)
   #spe=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=(1,2))
   #prof=M.sum(axis=0)
   #vmin=ma.min()
   #vmax=ma.max()
   #plt.clf() 
   #plt.subplot(2, 3, 1)
   #plt.imshow(ma,vmin=vmin,vmax=vmax)
   #plt.subplot(2, 3, 3)
   #plt.imshow(prof)
   #plt.subplot(2, 3, 2)

  

   # SHOW results (NOT ALGORITHMIC CODE)
   #spe2=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=(1,2))
   #plt.imshow(cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=0),vmin=vmin,vmax=vmax)
   #plt.subplot(2, 3, 6)
   #plt.imshow(syn.sum(axis=0))
   #plt.subplot(2, 3, 5)
   #plt.imshow(cube.data.sum(axis=0))
   #plt.gca().add_patch(plt.Rectangle((x_lb,y_lb),x_ub-x_lb+1,y_ub-y_lb+1,alpha=1, facecolor='none'))
   #plt.subplot(2, 3, 4)
   #plt.plot(spe,"b")
   #plt.plot(spe2,"r")
   #plt.show() 
   #plt.pause(0.01)

   # Return the clump parameters
   return res.x


def gauss_clumps_params():
   retval=dict()
   retval['threshold']=0.000001
   retval['few_deltas']=5
   retval['weight_deltas']=10
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
   params['beam_size']=abs(float(cube.meta['BMIN']))
   params['spe_res']=abs(float(cube.meta['CDELT3']))
   while not stop:
      theta=next_clump(cube,syn,params)
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
