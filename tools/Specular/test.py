import challa.workspace as ws
import challa.clumps as cl
import numpy as np
import timeit
from challa.cube import Cube

global M
global model

def f1():
   clumps._eval(M,model)

#ws.import_file("fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
ws.import_file("fits/Boom.cm.cln.fits")
elm=ws.elements()
#print elm
ndd=elm['Boom.cm.cln-0']
#ndd=elm['calibrated.ms.image.spectrum.J113740.6-010454.spw0.image-0']
params=cl.gc_default_params()
cb=Cube(ndd.data[0],ndd.meta)
theta=cl.gauss_clumps(cb,params)



