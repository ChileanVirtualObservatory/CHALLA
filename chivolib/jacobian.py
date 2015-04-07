import numpy as np
from sympy import *

# Script to calcule the jacobian of chi-square function


# Xv array with (x,y,v) values
# Yv array with intesity values
# w array with weights values
# Y gaussian function as symbolic expression
def chiJacobian():
        # Parameters which define a Gaussian clump
        veci=symbols('a b x0 y0 v0 phi sx sy sv dvx dvy')
        (a, b, x0, y0, v0, phi, sx, sy, sv, dvx, dvy) = veci

        # Variables of intensity Gaussian function
        (x, y, v) = symbols('x y v')
        (s0, sc, sa) = symbols('s_0 s_c s_a')
        (xmax1, xmax2, xmax3) = symbols('x_m y_m v_m')
	ymax = symbols('v_max')
	dBeam = symbols('\Delta_{beam}')
	dRes = symbols('\Delta_{res}')
        YFit = symbols('Y_Fit')
        YReal = symbols('Y_Real')
        w = symbols('w')
        #nf = symbols('nf')
        
        # Precision Matrix
        Lambda = np.array([
	   [cos(phi)**2/sx**2 + sin(phi)**2/sy**2 + dvx**2/sv**2, -sin(2*phi)/(2*sx**2) + sin(2*phi)/(2*sy**2) + dvx*dvy/sv**2, -dvx/sv**2],
	   [-sin(2*phi)/(2*sx**2) + sin(2*phi)/(2*sy**2) + dvx*dvy/sv**2, sin(phi)**2/sx**2 + cos(phi)**2/sy**2 + dvy**2/sv**2, -dvy/sv**2],
	   [-dvx/sv**2, -dvy/sv**2, 1/sv**2]])
        

        Lambda_phi=np.empty_like(Lambda)
        Lambda_sx=np.empty_like(Lambda)
        Lambda_sy=np.empty_like(Lambda)
        Lambda_sv=np.empty_like(Lambda)
        Lambda_dvx=np.empty_like(Lambda)
        Lambda_dvy=np.empty_like(Lambda)
        for i in range(3):
           for j in range(3):
              Lambda_phi[i][j]=diff(Lambda[i][j],phi)
              Lambda_sx[i][j]=diff(Lambda[i][j],sx)
              Lambda_sy[i][j]=diff(Lambda[i][j],sy)
              Lambda_sv[i][j]=diff(Lambda[i][j],sv)
              Lambda_dvx[i][j]=diff(Lambda[i][j],dvx)
              Lambda_dvy[i][j]=diff(Lambda[i][j],dvy)
        print "L_phi = " + str(Lambda_phi)
        print "L_sx = " + str(Lambda_sx)
        print "L_sy = " + str(Lambda_sy)
        print "L_sv = " + str(Lambda_sv)
        print "L_dvx = " + str(Lambda_dvx)
        print "L_dvy = " + str(Lambda_dvy)
        # x-mu vector
        v = np.array([x-x0, y-y0, v-v0])

        # Gaussian function
        Y = a * exp(-0.5*np.dot(v,np.dot(Lambda,v))) + b


	#s0 = 1
	#sc = 1
	#sa = 1
	#nf = Xv.size - 11
	# Jacobian of Y (parameters)
	Jy = np.array([diff(Y,a), diff(Y,b), diff(Y,x0), diff(Y,y0), diff(Y,v0), diff(Y,phi),
		diff(Y,sx), diff(Y,sy), diff(Y,sv), diff(Y,dvx), diff(Y,dvy)])
       
	chiJ = [0,0,0,0,0,0,0,0,0,0,0]

	for j in range(11):
        #        print "\\frac{\\partial{y}}{\\partial{"+str(veci[j])+"}}  & = "+str(latex(Jy[j])) + "\\\\"
		# Evaluates derivative respect to j-variable
		chiJ[j] += -2*w*(YReal-YFit)*Jy[j] + s0*exp(YFit-YReal)*Jy[j]

	#for j in range(11):
        #        print "\\frac{\\partial{\chi^2}}{\\partial{"+str(veci[j])+"}} &= "+str(latex(chiJ[j])) + "\\\\"
	# Sum extra terms for a and b derivatives
	(a, b) = symbols('a b')
	chiJ[0] += 2*sa*(a+b-ymax)
	chiJ[1] += 2*sa*(a+b-ymax)
        chiJ[2] += 2*sc*(xmax1-x)**2/dBeam**2
        chiJ[3] += 2*sc*(xmax2-y)**2/dBeam**2
        chiJ[4] += 2*sc*(xmax3-v)**2/dRes**2
	#print chiJ

chiJacobian()






	

