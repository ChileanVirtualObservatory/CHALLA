import numpy as np 
import sympy as sp
import pickle

"""
bbox function return the bounding box of an ellipsoid 
of the form (x-u)^T.A.(x-u)=k.
It returns (xinf,xsup,yinf,ysup,zinf,zsup)
"""
def lamb(func):
    return sp.lambdify((a,b,c,d,e,f,g,h,i), func, modules='numpy')

def bbox(A,k):
	###Load bounding functions
	fl=file("bbox.pickle","o")
	expr=pickle.load(fl)
	fl.close()

	###Transform Sympy expressions to Numpy evaluated functions
	functions=map(lamb,expr)

	###Applying functions
	args=A.reshape(1,9)
	(a,b,c,d,e,f,g,h,i)=tuple(args[0])
	box=list()
	for f in functions:
		box.append(f(a,b,c,d,e,f,g,h,i))

	###Ordering
	ret=np.min(box[0],box[1]),np.max(box[0],box[1]),
		np.min(box[2],box[3]),np.max(box[2],box[3]),
		np.min(box[4],box[5]),np.max(box[4],box[5]))
	return ret