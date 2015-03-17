import numpy as np
import copy

class Cube:
    """
    A spectral cube.
    """


    def __init__(self,data,meta):
        self.data=data
        self.meta=meta
        ra_value=float(meta['CRVAL1'])
        self.ra_delta=np.abs(float(meta['CDELT1']))
        ra_cpix=int(meta['CRPIX1']) -1
        ra_elms=int(meta['NAXIS1'])
        self.ra_axis=np.linspace(ra_value-ra_cpix*self.ra_delta,ra_value+(ra_elms-ra_cpix)*self.ra_delta, num=ra_elms)

        dec_value=float(meta['CRVAL2'])
        self.dec_delta=np.abs(float(meta['CDELT2']))
        dec_cpix=int(meta['CRPIX2']) -1
        dec_elms=int(meta['NAXIS2'])
        self.dec_axis=np.linspace(dec_value-dec_cpix*self.dec_delta,dec_value+(dec_elms-dec_cpix)*self.dec_delta, num=dec_elms)

        nu_value=float(meta['CRVAL3'])
        self.nu_delta=np.abs(float(meta['CDELT3']))
        nu_cpix=int(meta['CRPIX3']) -1   
        nu_elms=int(meta['NAXIS3'])
        self.nu_axis=np.linspace(nu_value-nu_cpix*self.nu_delta,nu_value+(nu_elms-nu_cpix)*self.nu_delta, num=nu_elms)

    def max(self):
        index=np.unravel_index(self.data.argmax(),self.data.shape)
        y=self.data[index]
        x=np.empty(3)
        #print index
        #print len(self.nu_axis)
        #print len(self.dec_axis)
        #print len(self.ra_axis)
        x[0]=self.nu_axis[index[0]]
        x[1]=self.dec_axis[index[1]]
        x[2]=self.ra_axis[index[2]]
        return (y,x)
     
#    def _check_dims(self,crval,crpix,naxis):
#        if np.abs(float(cb.meta['CDELT1'])) != np.abs(float(self.meta['CRDELT1'])) or np.abs(float(cb.meta['CDELT2'])) != np.abs(float(self.meta['CDELT2'])) or np.abs(float(cb.meta['CDELT3'])) != np.abs(float(self.meta['CDELT3'])):
#           print "WARNING: Different resolution substracting not implemented yet"
#           return
#        ra_val=float(cb.meta['CRVAL1'])
#        dec_val=float(cb.meta['CRVAL2'])
#        nu_val=float(cb.meta['CRVAL3'])
#        ra_ci=np.searchsorted(ra_axis, ra_val)
#        dec_ci=np.searchsorted(dec_axis, dec_val)
#        nu_ci=np.searchsorted(nu_axis, nu_val)
#        if ra_val!=ra_axis[ra_ci] or  dec_val!=dec_axis[dec_ci] or  nu_val!=nu_axis[nu_ci]:
#           print "WARNING: not aligned substracting not implemented yet"
#           return
#        if ra_ci + (int(cb.meta['NAXIS1']) - )
#        if naxis - crpix > olen - oci or crpix > oci:
#           print "WARNING: not aligned substracting not implemented yet"
#        return true

#    def substract(self,cb):
#        pos=self.find_position()

    def feature_space(self,center,window):
        print "subcube"
        #print center
        #print window
        ra_ci=np.argmin(np.abs(self.ra_axis-center[2]));
        ra_ui=np.argmin(np.abs(self.ra_axis-center[2]-window[2]));
        ra_li=np.argmin(np.abs(self.ra_axis-center[2]+window[2]));
        dec_ci=np.argmin(np.abs(self.dec_axis-center[1]));
        dec_ui=np.argmin(np.abs(self.dec_axis-center[1]-window[1]));
        dec_li=np.argmin(np.abs(self.dec_axis-center[1]+window[1]));
        nu_ci=np.argmin(np.abs(self.nu_axis-center[0]));
        nu_ui=np.argmin(np.abs(self.nu_axis-center[0]-window[0]));
        nu_li=np.argmin(np.abs(self.nu_axis-center[0]+window[0]));
        

        crval1=self.ra_axis[ra_ci]
        crval2=self.dec_axis[dec_ci]
        crval3=self.nu_axis[nu_ci]
        crpix1=ra_ci - ra_li 
        crpix2=dec_ci - dec_li 
        crpix3=nu_ci - nu_li 
        naxis1=ra_ui-ra_li +1 
        naxis2=dec_ui-dec_li +1
        naxis3=nu_ui-nu_li +1
        #print (crval3,crpix3,naxis3)
        #print ra_ui, ra_li
        ra_axis=np.linspace(crval1-crpix1*self.ra_delta,crval1+(naxis1-crpix1)*self.ra_delta, num=naxis1)
        dec_axis=np.linspace(crval2-crpix2*self.dec_delta,crval2+(naxis2-crpix2)*self.dec_delta, num=naxis2)
        nu_axis=np.linspace (crval2-crpix3*self.nu_delta,crval3+(naxis3-crpix3)*self.nu_delta, num=naxis3)

        adn=np.meshgrid(nu_axis,dec_axis,ra_axis, indexing='ij')
        X=np.empty((3,len(ra_axis)*len(dec_axis)*len(nu_axis)))
        X[0]=adn[0].ravel()
        X[1]=adn[1].ravel()
        X[2]=adn[2].ravel()
        print len(ra_axis)
        #y=self.data[nu_li:nu_ui,dec_li:dec_ui,ra_li:ra_ui].ravel()
        yidx=(nu_li,nu_ui,dec_li,dec_ui,ra_li,ra_ui)
        print "where"
        print yidx
        return X,yidx
        



#   def add(self,function,center,params,window):
       
