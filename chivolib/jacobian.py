import numpy as np
from sympy import *

# Script to calcule the jacobian of chi-square function


# Xv array with (x,y,v) values
# Yv array with intesity values
# w array with weights values
# Y gaussian function as symbolic expression
def chiJacobian(Xv,Yv,w,Y):
	s0 = 1
	sc = 1
	sa = 1
	nf = Xv.size - 11
	ymax = max(Yv)
	# Jacobian of Y (parameters)
	Jy = np.array([diff(Y,a), diff(Y,b), diff(Y,x0), diff(Y,y0), diff(Y,v0), diff(Y,phi),
		diff(Y,Dx), diff(Y,Dy), diff(Y,Dv), diff(Y,dvx), diff(Y,dvy)])
	chiJ = [0,0,0,0,0,0,0,0,0,0,0]

	for i in range(Xv.size):
		xv,yv,vv = Xv[i]
		Yfit = Y.subs(x,xv).subs(y,yv).subs(v,vv)
		for j in range(11):
			# Evaluates derivative respect to j-variable
			aux = Jy[j].subs(x,xv).subs.(y,yv).subs(v,vv) 
			chiJ[j] += -2*w[i]*(Y[i]-Yfit)*aux + s0*np.exp(Yfit-Y[i])*aux

	# Sum extra terms for a and b derivatives
	(a, b) = symbols('a b')
	chiJ[0] += 2*sa*(a+b-Ymax)
	chiJ[1] += 2*sa*(a+b-Ymax)
	return chiJ


# Parameters which define a Gaussian clump
(a, b, x0, y0, v0, phi, Dx, Dy, Dv, dvx, dvy) = symbols('a b x0 y0 v0 phi Dx Dy Dv dvx dvy')

# Variables of intensity Gaussian function
(x, y, v) = symbols('x y v')

# Precision Matrix
Lambda = np.array([
	[cos(phi)**2/Dx**2 + sin(phi)**2/Dy**2 + dvx**2/Dv**2, -sin(2*phi)/(2*Dx**2) + sin(2*phi)/(2*Dy**2) + dvx*dvy/Dv**2, -dvx/Dv**2],
	[-sin(2*phi)/(2*Dx**2) + sin(2*phi)/(2*Dy**2) + dvx*dvy/Dv**2, sin(phi)**2/Dx**2 + cos(phi)**2/Dy**2 + dvy**2/Dv**2, -dvy/Dv**2],
	[-dvx/Dv**2, -dvy/Dv**2, 1/Dv**2]])

# x-mu vector
v = np.array([x-x0, y-y0, v-v0])

# Gaussian function
Y = a * exp(-0.5*np.dot(v,np.dot(Lambda,v))) + b

# call chiJacobian() here...






	

