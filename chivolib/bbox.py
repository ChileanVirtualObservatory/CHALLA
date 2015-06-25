import numpy as np 
import sympy as sp
import pickle

"""
bbox function return the bounding box of an ellipsoid 
of the form (x-u)^T.A.(x-u)=k.
It returns (xinf,xsup,yinf,ysup,zinf,zsup)
"""
def lamb(func):
	a,b,c,d,e,f,g,h,i=sp.symbols("a b c d e f g h i")
	return sp.lambdify((a,b,c,d,e,f,g,h,i), func, modules='numpy')

if __name__=='bbox':
	###Load bounding expressions
	fl=file("bbox.pickle","r")
	expr=pickle.load(fl)
	fl.close()
	###Transform Sympy expressions to Numpy evaluated functions
	functions=map(lamb,expr)

def bbox(A,k=1):
	###Applying functions
	args=A.reshape(1,9)
	(a,b,c,d,e,f,g,h,i)=tuple(args[0])
	box=list()
	for j in range(6):
		box.append(functions[j](a,b,c,d,e,f,g,h,i))
	ret=(np.abs(box[0]),np.abs(box[2]),np.abs(box[4]))
	return ret


