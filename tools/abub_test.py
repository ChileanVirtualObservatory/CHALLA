import chivolib.workspace as ws
import chivolib.abub as bl
import numpy as np
import timeit
import cProfile
from chivolib.spectral import Cube

global M
global model

ws.import_file("fits/M100line.image.fits")
#ws.import_file("fits/combined-278000.fits")
#ws.import_file("fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
#ws.import_file("fits/calibrated.ms.line.spw0.source15.image.fits")
#ws.import_file("fits/Boom.cm.cln.fits")

elm=ws.elements()
#print elm
ndd=elm['M100line.image-0']
#ndd=elm['combined-278000-0']
#ndd=elm['calibrated.ms.image.spectrum.J113740.6-010454.spw0.image-0']
#ndd=elm['calibrated.ms.line.spw0.source15.image-0']
#ndd=elm['Boom.cm.cln-0']
#params=bl.gauss_clumps_params()
#params['threshold']=0.000
if ndd.meta['NAXIS']==4:
   cb=Cube(ndd.data[0],ndd.meta)
else:
   cb=Cube(ndd.data,ndd.meta)
#(vect,ener)=bl.bubble_fit(cb,1.5,0.5,plot=False)
(vect,ener)=bl.bubble_fit(cb,1.0,0.5)




