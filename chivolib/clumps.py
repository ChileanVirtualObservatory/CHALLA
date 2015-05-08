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
def chi2(model,features,values,w,value_max,feature_max,res_vect,s_vect):
   sys.stdout.write('.')
   sys.stdout.flush()
   su=values - gauss_eval(features,to_gauss(model))
   nf=len(su) - 11;
   t1 = (np.square(su)*w).sum()/nf
   t2 = (np.exp(-su)).sum()/nf
   t3 = np.square((model[2]-feature_max[0])/res_vect[0]) + np.square((model[3]-feature_max[1])/res_vect[1]) + np.square((model[4]-feature_max[2])/res_vect[2])
   t4 = np.square(model[0]+model[1] - value_max)
   #print "chi2"
   #print t1 + s_vect[0]*t2 + s_vect[1]*t3 + s_vect[2]*t4
   #print "mod"
   #print model
   return(t1 + s_vect[0]*t2 + s_vect[1]*t3 + s_vect[2]*t4)

def jac_chi2(model,features,values,w,value_max,feature_max,res_vect,s_vect):
   sys.stdout.write('*')
   sys.stdout.flush()

   # Unpack values
   (a, b, x0, y0, v0, phi, sx, sy, sv, dvx, dvy) = model
   (a,b,mu,L)=to_gauss(model)
   fit=gauss_eval(features,(a,b,mu,L))
   su=values - fit
   nf=len(su) - 11;
   basic_term = (-2*su*w + s_vect[0]*np.exp(-su))/nf
   #print "basic_term"
   #print basic_term
   exp_term = (fit - b)/a
   jaco=np.empty(11)
   # partial derivate w.r.t. a
   extra_sa = 2*s_vect[2]*(a + b - value_max)
   jaco[0]=(basic_term*exp_term).sum()  + extra_sa
   # partial derivate w.r.t. b
   jaco[1]=basic_term.sum() + extra_sa
   # partial derivate w.r.t. x0,y0 and v0
   C=np.empty_like(features)
   C[0]=features[0] - mu[0]
   C[1]=features[1] - mu[1]
   C[2]=features[2] - mu[2]
   V=-L.dot(C) # Derivate of the quadratic form
   extra_sc = 2*s_vect[1]*(mu - feature_max)/(res_vect*res_vect)
   # compose x0
   jaco[2]=(a*basic_term*exp_term*V[0]).sum() + extra_sc[0]
   # compose y0
   jaco[3]=(a*basic_term*exp_term*V[1]).sum() + extra_sc[1]
   # compose v0
   jaco[4]=(a*basic_term*exp_term*V[2]).sum() + extra_sc[2]
   ## partial derivate w.r.t. phi,sx,sy,sv,dvx,dvy
   sphi=np.sin(phi)
   s2phi=np.sin(2*phi)
   sphi2=np.square(sphi)
   cphi=np.cos(phi)
   c2phi=np.cos(2*phi)
   cphi2=np.square(cphi)
   sx2=np.square(sx)
   sy2=np.square(sy)
   sv2=np.square(sv)
   sx3=np.power(sx,3)
   sy3=np.power(sy,3)
   sv3=np.power(sv,3)
   dvx2=np.square(dvx)
   dvy2=np.square(dvy)
   D_phi=(1.0/sy2 - 1.0/sx2)*np.array([[2*sphi*cphi,c2phi,0],[c2phi,-2*sphi*cphi,0],[0,0,0]])
   V=-(C*(D_phi.dot(C))).sum(axis=0)/2.0
   jaco[5]=(a*basic_term*exp_term*V).sum()
   D_sx=(1.0/sx3)*np.array([[-2*cphi2,s2phi,0],[s2phi,-2*sphi2,0],[0,0,0]])
   V=-(C*(D_sx.dot(C))).sum(axis=0)/2.0
   jaco[6]=(a*basic_term*exp_term*V).sum()
   D_sy=(1.0/sy3)*np.array([[-2*sphi2,-s2phi,0],[-s2phi,-2*cphi2,0],[0,0,0]])
   V=-(C*(D_sy.dot(C))).sum(axis=0)/2.0
   jaco[7]=(a*basic_term*exp_term*V).sum()
   D_sv=(2.0/sv3)*np.array([[-dvx2,-dvx*dvy,dvx],[-dvx*dvy,-dvy2,dvy],[dvx,dvy,-1]])
   V=-(C*(D_sv.dot(C))).sum(axis=0)/2.0
   jaco[8]=(a*basic_term*exp_term*V).sum()
   D_dvx=(1.0/sv2)*np.array([[2*dvx,dvy,-1],[dvy,0,0],[-1,0,0]])
   V=-(C*(D_dvx.dot(C))).sum(axis=0)/2.0
   jaco[9]=(a*basic_term*exp_term*V).sum()
   D_dvy=(1.0/sv2)*np.array([[0,dvx,0],[dvx,2*dvy,-1],[0,-1,0]])
   V=-(C*(D_dvy.dot(C))).sum(axis=0)/2.0
   jaco[10]=(a*basic_term*exp_term*V).sum()
   #print "jaco"
   #print jaco
   return jaco
 

def standarize((features,values,w,value_max,feature_max,res_vect,s_vect)):
   s_features=np.empty_like(features)
   s_xmax=np.abs(features[0]).max()
   s_ymax=np.abs(features[1]).max()
   s_zmax=np.abs(features[2]).max()
   s_features[0]=features[0]/s_xmax
   s_features[1]=features[1]/s_ymax
   s_features[2]=features[2]/s_zmax
   s_features_max=np.empty_like(feature_max)
   s_features_max[0]=feature_max[0]/s_xmax
   s_features_max[1]=feature_max[1]/s_ymax
   s_features_max[2]=feature_max[2]/s_zmax
   s_values=values/value_max
   s_res_vect=np.empty_like(res_vect)
   s_res_vect[0]=res_vect[0]/s_xmax
   s_res_vect[1]=res_vect[1]/s_ymax
   s_res_vect[2]=res_vect[2]/s_zmax
   args=(s_features,s_values,w,1,s_features_max,s_res_vect,s_vect)
   #(a,b,x0,y0,v0,phi,sx,sy,sv,dvx,dvy)
   tr=np.array([value_max,value_max,s_xmax,s_ymax,s_zmax,1,s_xmax,s_ymax,s_zmax,s_zmax/s_xmax,s_zmax/s_ymax])
   return (args,tr)

def next_clump(cube,syn,params):
   # Non-blocking plot 
   plt.ion() 
   plt.clf() 
   
   (value_max,feature_max) = cube.max()
   b=params['rms'] # need to be a local minima... not good
   a=value_max - b
   
   # Initial guess: position and orientation
   x0=feature_max[0]
   y0=feature_max[1]
   v0=feature_max[2]
   phi=0
   
   # Initial guess: variances
   res_vect=np.array([params['beam_size'],params['beam_size'],params['spe_res']])
   s_vect=np.array([params['s0'],params['sc'],params['sa']])
 
   # Initial guess: redshift
   dvalp=0
   dvdel=0
  
   # Compute the weight vector and the feature space
   w_sigmas=params['weight_deltas']*res_vect # several times the resolution
   (features,sc_index) = cube.feature_space(feature_max,2*w_sigmas) 
   w_shape=(1,0,x0,y0,v0,0,w_sigmas[0],w_sigmas[1],w_sigmas[2],0,0)
   w=gauss_eval(features,to_gauss(w_shape),False)
 
   #Plot current cube
   plt.subplot(2, 3, 4)
   plt.imshow(cube.stack())
   rect=plt.Rectangle((sc_index[0],sc_index[2]),sc_index[1]-sc_index[0]+1,sc_index[3]-sc_index[2]+1,alpha=1, facecolor='none')
   plt.gca().add_patch(rect)
   
   #Plot current subcube
   plt.subplot(2, 3, 1)
   plt.imshow(cube.stack(sc_index))

   
   # Ravel the values
   values=cube.ravel(sc_index)
   
   idx=(values-a/2).argmin()
   sq2l2=np.sqrt(2*np.log(2))
   sx=np.abs(features[0][idx]-feature_max[0])/sq2l2
   sy=np.abs(features[1][idx]-feature_max[1])/sq2l2
   sv=np.abs(features[2][idx]-feature_max[2])/sq2l2
   # Compile first guess
   guess=np.array([a,b,x0,y0,v0,phi,sx,sy,sv,dvalp,dvdel])
   print "GUESS"
   print guess
   
   # Pack all args of chi2
   chi2_args=(features,values,w,value_max,feature_max,res_vect,s_vect)
   

   print "grad error", check_grad(chi2,jac_chi2,guess,features,values,w,value_max,feature_max,res_vect,s_vect)
   #print "PREV= ", chi2_args
   # Standarize everything
   #(std_args,tr_vect)=standarize(chi2_args)
   #print "STAND= ", std_args
   #print "tr= ", tr_vect
   #s_guess=guess/tr_vect
   # OPTIMIZE
   #res = minimize(chi2,guess,jac=jac_chi2,method='CG',args=chi2_args)
   res=fmin_bfgs(chi2, guess,fprime=jac_chi2, args=chi2_args)
   #res=fmin_bfgs(chi2, guess,args=chi2_args)
   #res = minimize(chi2,guess,jac=jac_chi2,method='BFGS',args=chi2_args,tol=1e-30)
   print
   print "res =", res
   print
   print "clump =", res  #*tr_vect
   sys.stdout.flush()
   # Clump values
   #print "AND NOW THE GAUSS"
   #print to_gauss(res.x) 
   val_fit=gauss_eval(features,to_gauss(res),False)
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
   plt.show() 
   plt.pause(0.01)

   # Return the clump parameters
   return res


def gauss_clumps_params():
   retval=dict()
   retval['threshold']=0.000001
   retval['weight_deltas']=20
   retval['s0']=1.0
   retval['sc']=1.0
   retval['sa']=1.0
   #retval['rms']=0.00002
   return retval

def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)

def gauss_clumps(orig_cube,params):
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
   print params['rms']
   while not stop:
      theta=next_clump(cube,syn,params)
      C.append(theta)
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
   print C
   return C
