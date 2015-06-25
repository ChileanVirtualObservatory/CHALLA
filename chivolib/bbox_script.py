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
of functions that gerate the box limits.
"""

def lamb(func):
    return sp.lambdify((a,b,c,d,e,f,g,h,i), func, modules='numpy')

(a,b,c,d,e,f,g,h,i)=sp.symbols('a b c d e f g h i')
(x,y,z)=sp.symbols('x y z')
A=np.array([[a,b,c],[d,e,f],[g,h,i]])
X=np.array([[x],[y],[z]])
quad=np.dot(np.dot(X.T,A),X)[0,0]
eq=quad-1
eqx=sp.diff(eq,x)
eqy=sp.diff(eq,y)
eqz=sp.diff(eq,z)

###Solve for all bounds ###
solx=sp.solve([eq,eqy,eqz],[x,y,z])
soly=sp.solve([eq,eqx,eqz],[x,y,z])
solz=sp.solve([eq,eqx,eqy],[x,y,z])

###Simplify the bounds ###
x1=sp.simplify(solx[0][0])
x2=sp.simplify(solx[1][0])
y1=sp.simplify(soly[0][1])
y2=sp.simplify(soly[1][1])
z1=sp.simplify(solz[0][2])
z2=sp.simplify(solz[1][2])

###Convert symbolic expressions to functions evaluated with numpy###
x1,x2,y1,y2,z1,z2=map(lamb,(x1,x2,y1,y2,z1,z2))
functions=(x1,x2,y1,y2,z1,z2)

###Serializing...###
f=file("bbox.pickle","w")
pickle.dump(functions,f)
f.close()

