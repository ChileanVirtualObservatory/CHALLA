import challa.workspace as ws
import challa.clumps as clumps
import numpy as np

ws.import_file("../fits/logfile_alma_hatlas_cycle1_inc-z_beye.fits")
elm=ws.elements()
ndd=elm['logfile_alma_hatlas_cycle1_inc-z_beye-1']
model=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1)
val=clumps._eval(np.random.rand(3,100),model)
print val
print len(val)

