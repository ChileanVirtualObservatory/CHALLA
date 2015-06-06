import numpy as np
from collections import deque
from scipy.optimize import fmin_bfgs,check_grad
import copy
import matplotlib.pyplot as plt
import sys
from cube import *

class GaussClumps:

   def __init__(self):
      # Set the very 
      self.defaultParams()
   
   def defaultParams(self):
      self.par=dict()
      # Spectral Resoluion in pixels (smoothing function)
      self.par['VELORES']=2.0
      # Beam resoluion in pixels (smoothing function)
      self.par['FWHMBEAM']=2.0
      # The maximum allowed number of failed fits between succesful fits.
      self.par['MAXSKIP']=10
      # Maximum Clumps
      self.par['MAXCLUMPS']=sys.maxint
      # The iterative process ends when "npad" consecutive clumps all had peak
      # values below "peak_thresh" or all had areas below "area_thresh".
      self.par['NPAD']=10
      # The lower threshold for clump peaks to a user-specified multiple of the RMS noise.
      self.par['THRESH']=2.0
      # The lower threshold for clump area to a user-specified number of pixels.
      self.par['MINPIX']=3
      # The lowest value (normalised to the RMS noise level) at which
      # model Gaussians should be evaluated. 
      self.par['MODELMIN']=0.5
      # The max allowed fraction of bad pixels in a clump. 
      self.par['MAXBAD']=0.05
      # No.of standard deviations at which to reject peaks
      self.par['NSIGMA']=3.0
      # But reject peaks only if at least NPEAKS were found 
      self.par['NPEAKS']=9
   
   def fit(self,cube,rms=-1.0,verbose=False):
      # Set the RMS, or automatically find an estimate for it
      if rms < 0.0:
         rms=cube.estimate_rms()
      self.par['RMS']=rms
      
      # Unpack used parameters
      npeaks=self.par['NPEAKS']
      mlim=self.par['MODELMIN']
      peak_thresh=self.par['THRESH']
      area_thresh=self.par['MINPIX']
      maxclump=self.par['MAXCLUMPS']
      npad=self.par['NPAD']

      # Copy the supplied cube into a work cube which will hold the
      # residuals remaining after subtraction of the fitted Gaussians. 
      self.res=copy.deepcopy(cube)

      # Initialise the number of clumps found so far.
      iclump = 0

      # Indicate that no peaks have been found below the lower threshold for clump
      # peak values, or below the lower area threshold. 
      peaks_below = 0
      area_below = 0
 
      # Initialise the variables used to keep track of the mean and standard
      # deviation of the most recent "npeak" fitted peak values. 
      mean_peak = 0.0;
      sigma_peak = 0.0; 
      # The value most recently added to "peaks"
      new_peak = 0.0;
      # Sum of the values in "peaks"
      sum_peak = 0.0
      # Sum of the squares of the values in "peaks" 
      sum_peak2 = 0.0

      # Number of pixels contributing to the clump
      area=0

      # Iterations performed so far 
      niter = 0
      iterate = True
      # No. of failed fits since last good fit 
      nskip = 0
      # Sum of the values in all the used clumps so far 
      sumclumps = 0.0
      # Sum of the supplied data values 
      sumdata = np.nan
      
      # peaks contains the last npeaks... 
      peaks=np.zeros(npeaks)
      
      # Loop round fitting a gaussian to the largest remaining peak in the
      # residuals array. */
      while iterate:
         # Report the iteration number to the user if required.
         niter+=1
         if verbose:
            print "Iteration:", niter
         # Find the cube index of the element with the largest value in the residuals cube.
         # imax: Index of element with largest residual
         (allbad,imax,sumdata) = self.findMax(sumdata)
         # imax: Index of element with largest residual
 
         # Finish iterating if all the residuals are bad, or if too many iterations
         # have been performed since the last succesfully fitted clump. 
         if allbad:
            iterate = False;
            niter-=1
            if verbose:
               print "There are no good pixels left to be fitted."
               continue
         elif  nskip > maxskip:
            iterate = False;
            niter-=1
            if verbose:
               print "The previous",maxskip,"fits were unusable."
               continue
         # If not, make an initial guess at the Gaussian clump parameters centred on the current peak.
         guess=self.setInit(imax,niter)
         
         # Find the best fitting parameters, starting from the above initial guess.
         (found,clump,chisq)=self.optimize(imax,guess)
         # If no fit could be performed, then found = False
         if found:
            # Skip this fit if we have an estimate of the standard deviation of the
            # "npeaks" most recent clump peak values, and the peak value of the clump
            # just fitted is a long way (more than NSIGMA standard deviations) from the
            # peak value of the previously fitted clump. Also skip it if the peak
            # value is less than the "mlim" value.
            if (peaks.size == 0 or iclump < npeak or np.abs(clump[0] - new_peak) < nsig*sigma_peak ) and clump[0] > mlim]:

               # Record the new peak value for use with the next peak, and update the
               # standard deviation of the "npeak" most recent peaks. These values are
               # stored cyclically in the "peaks" array. */
               if peaks.size > 0:
                  np.roll(peaks,1)
                  new_peak = clump[0]
                  old_peak = peaks[0]
                  peaks[0]=new_peak
                  sum_peak += new_peak - old_peak
                  sum_peak2 += new_peak*new_peak - old_peak*old_peak
                  if sum_peak2 < 0.0:
                     sum_peak2 = 0.0
                  mean_peak = sum_peak/npeaks
                  sigma_peak = np.sqrt(sum_peak2/npeaks - mean_peak*mean_peak)
               
               # Increment the number of peaks found. 
               iclump+=1

               # Reset the number of failed fits since the last good fit. */
               nskip = 0;

               # Remove the model fit (excluding the background) from the residuals.
               # This also creates data values asociated with the clumps
               # The standard deviation of the new residuals is returned. */
               (result,area,sumclumps)=self.updateResults(clump,imax,mean_peak,sumclumps) # check check check
               #cupidGCUpdateArrays( type, res, ipd, el, ndim, dims,
               #                    x, rms, mlim, imax, peak_thresh, slbnd,
               #                    &ret, iclump, excols, mean_peak,
               #                    maxbad, &area, &sumclumps, status );

               # TODO: implement this!
               # Display the clump parameters on the screen if required. */
               #cupidGCListClump( iclump, ndim, x, chisq, slbnd, rms, status );

               # If this clump has a peak value which is below the threshold, increment
               # the count of consecutive clumps with peak value below the threshold.
               # Otherwise, reset this count to zero.
               if clump[0] < peak_thresh:
                  peaks_below+=1
               else:
                  peaks_below=0

               # If this clump has an area which is below the threshold, increment
               # the count of consecutive clumps with area below the threshold.
               # Otherwise, reset this count to zero. 
               if area < area_thresh:
                  area_below+=1
               else:
                  area_below=0

               # If the maximum number of clumps have now been found, exit.*/
               if iclump == maxclump:
                  iterate = False
                  if verbose:
                     print "The specified maximum number of clumps (",maxclump,") have been found."

               # If the integrated data sum in the fitted gaussians exceeds or equals
               # the integrated data sum in the input, exit. 
               elif sumclumps >= sumdata:
                  iterate = False
                  if verbose:
                     print "The total data sum of the fitted Gaussians (",sumclumps") has reached the total data sum in the supplied data (",sumdata,")."
 
               # If the count of consecutive peaks below the threshold has reached
               # "Npad", terminate.
               elif peaks_below == npad:
                  iterate = False
                  if verbose:
                     print "The previous",npad,"clumps all had peak values below the threshold."

               # If the count of consecutive clumps with area below the threshold has reached
               # "Npad", terminate.
               elif area_below == npad:
                  iterate = False
                  if verbose:
                     print "The previous ",npad," clumps all had areas below the threshold."
               
            # If the peak value fitted is very different from the previous fitted peak
            # value, set the residuals array element bad in order to prevent the
            # algorithm from trying to fit a peak to the same pixel again. 
            else:
               res[imax]=np.nan;
               new_peak = 0.5*(new_peak + clump[0]);
               nskip++;
               if verbose:
                  print "Clump rejected due to aberrant peak value. Ignoring Pixel..."

         # Tell the user if no clump could be fitted around the current peak
         # pixel value 
         else:
            nskip++;
            if verbose:
               print "No clump fitted (optimization falied). Ignoring Pixel..."
            # Set the specified element of the residuals array bad if no fit was
            # performed. This prevents the any subsequent attempt to fit a Gaussian
            # to the same peak value.
            res[imax]=np.nan;
      if verbose:
         print "GaussClump finished normally"
         # TODO: Usable Clumps
         ## Tell the problems with the clumps. */
         #if nclump == 0:
         #   print "No usable clumps found."
         #if iclump - nclump >= 1:
         #   print iclump - nclump,"clump(s) rejected because they touch an edge of the data array."
         ## Tell the user how many iterations have been performed (i.e. how many
         ## attempts there have been to fit a Gaussian peak
         #if niter == 1:
         #   print "No fit attempted."
         #else:
         #   print "Fits attempted for ",iclump," candidate clumps (",niter-iclump," failed)."
      return result

       

def chi2(model,features,values,w,value_max,feature_max,res_vect,s_vect):
   sys.stdout.write('.')
   sys.stdout.flush()
   su=values - gauss_eval(features,to_gauss(model))
   nf=len(su) - 11;
   t1 = (np.square(su)*w).sum()/nf
   t2 = (np.exp(-su)).sum()/nf
   t3 = np.square((model[2]-feature_max[0])/res_vect[0]) + np.square((model[3]-feature_max[1])/res_vect[1]) + np.square((model[4]-feature_max[2])/res_vect[2])
   t4 = np.square(model[0]+model[1] - value_max)
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
   # Compute basic term
   basic_term = (-2*su*w + s_vect[0]*np.exp(-su))/nf
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
   return jaco

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
