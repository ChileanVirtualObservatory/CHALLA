import math
import numpy as np
# Model: a,b,x0,y0,v0,dvx,dvy,Dx,Dy,Dv,phi

def _find_max(image):
   num=image.max()
   index=np.unravel_index(image.argmax(),image.shape)
   return (num,index)

# Older and Slower implementation (like in Stutzki)
#def _eval(M,(a,b,alp0,del0,v0,dvalp,dvdel,Dalp,Ddel,Dv,phi)):
#   sphi=math.sin(phi)
#   cphi=math.cos(phi)
#   s2phi=math.sin(2*phi)
#   Ga=math.pow(cphi/Dalp,2) + math.pow(sphi/Ddel,2)
#   Gc=math.pow(sphi/Dalp,2) + math.pow(cphi/Ddel,2)
#   Gb=-s2phi/(2*np.power(Dalp,2)) + s2phi/(2*np.power(Ddel,2)) 
#   (X,Y,V)=(M[0] - alp0,M[1] - del0,M[2] - v0)
#   quad=Ga*np.power(X,2) + 2*Gb*X*Y + Gc*np.power(Y,2)
#   shift=np.power((V - X*dvalp - Y*dvdel)/Dv,2)
#   retval=a*np.exp(-quad/2 - shift/2) + b
#   return retval


def _model2gauss(((a,b,alp0,del0,v0,phi,Dalp,Ddel,Dv,dvalp,dvdel)))
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
   return (a,b,

def _eval(X,M):
   sphi2=math.pow(math.sin(phi),2)
   cphi2=math.pow(math.cos(phi),2)
   s2phi=math.sin(2*phi)
   Dalp2=math.pow(Dalp,2)
   Ddel2=math.pow(Ddel,2)
   Dv2=math.pow(Dv,2)
   Ga=cphi2/Dalp2 + sphi2/Ddel2 + math.pow(dvalp,2)/Dv2
   Gd=sphi2/Dalp2 + cphi2/Ddel2 + math.pow(dvdel,2)/Dv2
   Gb=-s2phi/(2*Dalp2) + s2phi/(2*Ddel2) + dvalp*dvdel/Dv2
   Gc=-dvalp/Dv2
   Ge=-dvdel/Dv2
   Gf=1.0/Dv2
   G=np.array([[Ga,Gb,Gc],[Gb,Gd,Ge],[Gc,Ge,Gf]])
   C=np.empty_like(M)
   C[0]=X[0] - alp0
   C[1]=X[1] - del0
   C[2]=X[2] - v0
   V=C*(G.dot(C))
   quad=V.sum(axis=0)
   retval=a*np.exp(-quad/2) + b
   return retval



#def chi(alpha


#def gauss_clumps(ndd,sigma):
#   nmax,imax=_find_max(ndd.data)
