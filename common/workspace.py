
from astropy.io import fits
from astropy.nddata import NDData
import numpy as np
import os
from astropy import log



def ws_new(name):
   ws=dict()
   ws['name']=name


def ws_fits_consumer(ws):

def ws_hdf5_consumer(ws):
   log.setLevel('INFO')
def ws_votable_consumer(ws):
def ws_ascii_consumer(ws):
   


def ws_import(ws,path):
   filename=os.path.basename(path)
   name,ext=os.path.splitext(filename)
   if ext = '.fits':
      _ws_fits_consumer(ws,path,name)
   elif ext = '.hdf5':
      _ws_hdf5_consumer(ws,path,name)
   elif ext = '.xml':
      _ws_votable_consumer(ws,path,name)
   else:
      _ws_ascii_consumer(ws,path,name)
   
