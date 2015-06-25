import numpy as np 
import sympy as sp
import pickle

"""
Script to generate a general expression to the bounding box
of an ellipsoid of the form (x-u)^T.A.(x-u)=1.
Because we only want to know the dimensions of the bounding
box, i.e, displacement from center point, we center the ellipsoid
at u=(0,0,0) for simplicity.
It generates a file called "bbox.pickle" that serialize a tuple
of expressions for the box limits.
"""

(a,b,c,d,e,f,g,h,i)=sp.symbols('a b c d e f g h i')
(x,y,z)=sp.symbols('x y z')
A=np.array([[a,b,c],[d,e,f],[g,h,i]])
X=np.array([[x],[y],[z]])
quad=np.dot(np.dot(X.T,A),X)[0,0]
eq=quad-1
eqx=sp.diff(eq,x)
eqy=sp.diff(eq,y)
eqz=sp.diff(eq,z)

###Solve for all bounds 
solx=sp.solve([eq,eqy,eqz],[x,y,z])
soly=sp.solve([eq,eqx,eqz],[x,y,z])
solz=sp.solve([eq,eqx,eqy],[x,y,z])

###Simplify the bounds 
xlim=sp.simplify(solx[0][0])
ylim=sp.simplify(soly[0][1])
zlim=sp.simplify(solz[0][2])
obj=(xlim,ylim,zlim)

###Serializing...
f=file("bbox.pickle","w")
pickle.dump(obj,f)
f.close()



