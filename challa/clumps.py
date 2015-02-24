import math
import numpy as np
# Model: a,b,x0,y0,v0,dvx,dvy,Dx,Dy,Dv,phi

def _find_max(image):
   num=image.max()
   index=np.unravel_index(image.argmax(),image.shape)
   return (num,index)

def _eval(M,(a,b,alp0,del0,v0,dvalp,dvdel,Dalp,Ddel,Dv,phi)):
   sphi=math.sin(phi)
   cphi=math.cos(phi)
   s2phi=math.sin(2*phi)
   #TODO: Convert Dx and Dy to sigmas (they are FWHM i think...)
   Ga=(math.pow(cphi/Dalp,2) + math.pow(sphi/Ddel,2))/2.0
   Gc=(math.pow(sphi/Dalp,2) + math.pow(cphi/Ddel,2))/2.0
   Gb=-s2phi/np.power(2*Dalp,2) + s2phi/np.power(2*Ddel,2)
   #G=np.array([[Ga, Gb],[Gb, Gc]])
   #grad=np.array([dvalp,dvdel])
   (X,Y,V)=(M[0] - alp0,M[1] - del0,M[2] - v0)
   quad=Ga*np.power(X,2) + 2*Gb*X*Y + Gc*np.power(Y,2)
   shift=4*math.log(2)*np.power((V - X*dvalp - Y*dvdel)/Dv,2)
   retval=a*np.exp(-quad - shift) + b
   return retval

#def chi(alpha


#def gauss_clumps(ndd,sigma):
#   nmax,imax=_find_max(ndd.data)
