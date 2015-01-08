import challa.workspace as ws

ws.import_file("../fits/logfile_alma_hatlas_cycle1_inc-z_beye.fits")
for elm in ws.elements().keys():
   print elm
   ws.send(elm)
#ws.send('p33SO2-obs1-2')

