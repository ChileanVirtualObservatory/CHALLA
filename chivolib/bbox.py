import numpy as np 
import sympy as sp
import pickle

"""
bbox funtion return the bounding box of an ellipsoid 
of the form (x-u)^T.A.(x-u)=k.
"""

def bbox(A,k):
	f=file("bbox.pickle","o")
	functions=pickle.load(f)
	f.close()
