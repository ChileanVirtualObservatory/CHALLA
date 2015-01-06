
from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table
import numpy as np
import os
from astropy import log
from astropy.vo.samp import SAMPIntegratedClient
from urlparse import urlparse

def create(name):
   ws=dict()
   ws['workspace']=name
   return ws

def send(ws,name):
   _send(ws[name],name)

#This function was copied from a astropy development branch, please remove it when highlevel commands become available
def _send(data, name, destination='all', timeout=10, hub=None):
    """
    Send data to SAMP clients.

    Parameters
    ----------
    data : `~astropy.table.table.Table` or `~astropy.nddata.nddata.NDData` or `~numpy.ndarray` or `~astropy.io.fits.PrimaryHDU` or `~astropy.io.fits.ImageHDU` or `~astropy.io.fits.BinTableHDU` or `~astropy.io.fits.TableHDU`
        The data to send over SAMP
    name : str, optional
        The name of the dataset to use in other SAMP clients
    destination : str, optional
        The client to send the data to. By default, the data is broadcast to
        all SAMP clients. You can find the full list of available clients, use
        the :func:`~astropy.vo.samp.high_level.list_clients` function. As a
        convenience, you can also use ``'ds9'``, ``'topcat'``, and ``aladin'``
        and :func:`~astropy.vo.samp.high_level.send` will try and identify the
        correct client.
    timeout : int, optional
        The timeout for the request.
    hub : `~astropy.vo.samp.hub.SAMPHubServer`, optional
        The hub to send the data through
    """

    message = {}
    output_file = tempfile.NamedTemporaryFile()

    if isinstance(data, Table):

        data.write(output_file, format='votable')
        message['samp.mtype'] = "table.load.votable"

    elif isinstance(data, NDData):

        data.write(output_file, format='fits')
        message['samp.mtype'] = "image.load.fits"

    elif isinstance(data, np.ndarray):

        if data.dtype.fields is None:
            fits.writeto(output_file, data)
            message['samp.mtype'] = "image.load.fits"
        else:
            data = Table(data)
            data.write(output_file, format='votable')
            message['samp.mtype'] = "table.load.votable"

    elif isinstance(data, (fits.ImageHDU, fits.PrimaryHDU)):

        data.writeto(output_file)
        message['samp.mtype'] = "image.load.fits"

    elif isinstance(data, (fits.BinTableHDU, fits.TableHDU)):

        data.writeto(output_file)
        message['samp.mtype'] = "table.load.fits"

    else:

        raise TypeError("Unrecognized data type: {0}".format(type(data)))

    message['samp.params'] = {"url": urljoin('file:', output_file.name),
                              "name": name}

    client = SAMPIntegratedClient()
    client.connect(hub=hub)
    declare_metadata(client)

    if destination == 'all':
        for c in client.get_subscribed_clients(message['samp.mtype']):
            client.call_and_wait(c, message, timeout=str(timeout))
    elif destination in ['ds9', 'topcat', 'aladin']:
        clients = list_clients()
        for target_client in clients:
            name = target_client['name']
            if destination in name.lower():
                client.call_and_wait(str(target_client['id']), message,
                                     timeout=str(timeout))
    else:
        client.call_and_wait(destination, message, timeout=str(timeout))

    client.disconnect()

    output_file.close()

def _fits_consumer(ws,path,name):
#TODO: Support more filetypes
   wname=ws['workspace']
   log.info("Loading "+name+".fits into "+wname)
   hdulist = fits.open(path)
   counter=0
   for hdu in hdulist:
      if isinstance(hdu,fits.PrimaryHDU) or isinstance(hdu,fits.ImageHDU):
         log.info("Processing HDU "+str(counter)+" (image)")
         #TODO: check for WCS data...
         ndd=NDData(hdu.data,meta=hdu.header)
         ide=name+"-"+str(counter)
         ws[ide]=ndd
      elif isinstance(hdu,fits.BinTableHDU) or isinstance(hdu,fits.TableHDU):
         log.info("Processing HDU "+str(counter)+" (Table)")
          #TODO: check for WCS data...
         ntt=Table(hdu.data,meta=hdu.header)
         ide=name+"-"+str(counter)
         ws[ide]=ntt
      else:
         log.warning("HDU type not recognized, ignoring "+hdu.name+" ("+counter+")")
      counter+=1

def _hdf5_consumer(ws,path,name):
   log.warning("HDF5 format not supported yet. Ignoring file "+name+".hdf5")
def _votable_consumer(ws,path,name):
   log.warning("VOTable format not supported yet. Ignoring file "+name+".xml")
def _ascii_consumer(ws,path,name):
   log.warning("ASCII format not supported yet. Ignoring file "+name)


def file_import(ws,path):
   filename=os.path.basename(path)
   name,ext=os.path.splitext(filename)
   if ext == '.fits':
      _fits_consumer(ws,path,name)
   elif ext == '.hdf5':
      _hdf5_consumer(ws,path,name)
   elif ext == '.xml':
      _votable_consumer(ws,path,name)
   else:
      _ascii_consumer(ws,path,name)
   
