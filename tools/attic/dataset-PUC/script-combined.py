from challa import vu
import random
import math
import sys

rvel=150.0
temp=300.0
molist=('COv=0','13COv=0','C18O','C17O','13C18O','NH2','N2H+v=0','CNv=0','HCNv=0','HNCv=0','H2CN','CSv=0','CCS','H2S','H2CS','SO2v=0','H2CO','HCO+v=0','HC3Nv=0','HC5Nv=0','CH3OHvt=0')

log=open('script-combined.log', 'w')
dbpath="../ASYDO"

univ=vu.Universe(log)
for mol in molist:
   univ.create_source('combined-'+mol,0.0,0.0)
   s_x=random.uniform(0.2, 1)
   s_y=random.uniform(0.2, 1)
   s_f=random.uniform(2000, 8000)
   z_x=random.uniform(-10, 10)
   z_y=random.uniform(-10, 10)
   #curt=random.uniform(-3, 3)
   angle=random.uniform(0,math.pi)
   model=vu.IMCM(log,dbpath,mol,temp,(s_x,s_y,s_f,angle),(z_x,z_y))
   model.set_radial_velocity(rvel)
   univ.add_component('combined-'+mol,model)
for i in range(49):
   fcenter=(276 + 2*i)*1000 
    name, pos,res,pix,crpix
   cube=univ.gen_cube('combined',(0.0,0.0,fcenter))
   univ.save_cube(cube,'combined-'+str(fcenter)+'.fits')

