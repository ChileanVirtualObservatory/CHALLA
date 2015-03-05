

class Cube:
    """
    A spectral cube.
    """


    def __init__(self,data,meta):
        self.data=data
        self.meta=meta
        ra_value=float(meta['CRVAL1'])
        ra_delta=math.abs(float(meta['CDELT1']))
        ra_cpix=float(meta['CRPIX1']) -1
        ra_elms=float(meta['NAXIS1'])
        self.ra_axis=np.linspace(np.linspace(ra_value-ra_cpix*ra_delta,ra_value+(ra_elms-ra_cpix)*ra_delta, num=ra_elms))

        dec_value=float(meta['CRVAL2'])
        dec_delta=math.abs(float(meta['CDELT2']))
        dec_cpix=float(meta['CRPIX2']) -1
        dec_elms=float(meta['NAXIS2'])
        self.dec_axis=np.linspace(np.linspace(dec_value-dec_cpix*dec_delta,dec_value+(dec_elms-dec_cpix)*dec_delta, num=dec_elms))

        nu_value=float(meta['CRVAL3'])
        nu_delta=math.abs(float(meta['CDELT3']))
        nu_cpix=float(meta['CRPIX3']) -1   
        nu_elms=float(meta['NAXIS3'])
        nu_axis=np.linspace(np.linspace(nu_value-nu_cpix*nu_delta,nu_value+(nu_elms-nu_cpix)*nu_delta, num=nu_elms))

    def max(self):
        index=np.unravel_index(image.argmax(),image.shape)
        y=self.data[index]
        x=np.empty(3)
        x[0]=self.ra_axis[index[0]]
        x[1]=self.dec_axis[index[1]]
        x[2]=self.nu_axis[index[2]]
        return (y,x)
        

    def subcube(self,center,window):
        ra_ci=np.argmin(abs(self.ra_axis-center[0]));
        ra_ui=np.argmin(abs(self.ra_axis-center[0]+window[0]));
        ra_li=np.argmin(abs(self.ra_axis-center[0]-window[0]));
        dec_ci=np.argmin(abs(self.dec_axis-center[1]));
        dec_ui=np.argmin(abs(self.dec_axis-center[1]+window[1]));
        dec_li=np.argmin(abs(self.dec_axis-center[1]-window[1]));
        nu_ci=np.argmin(abs(self.nu_axis-center[2]));
        nu_ui=np.argmin(abs(self.nu_axis-center[2]+window[2]));
        nu_li=np.argmin(abs(self.nu_axis-center[2]-window[2]));
        new_data=self.data[ra_li:ra_ui,dec_li:dec_ui,nu_li:nu_li]
        new_meta=copy.deepcopy(self.meta)
        new_meta['CRVAL1']=self.ra_axis[ra_ci]
        new_meta['CRVAL2']=self.dec_axis[dec_ci]
        new_meta['CRVAL3']=self.nu_axis[nu_ci]
        new_meta['CRPIX1']=ra_ci - ra_li + 1
        new_meta['CRPIX1']=dec_ci - dec_li + 1
        new_meta['CRPIX1']=nu_ci - nu_li + 1
        new_meta['NAXIS1']=new_data.shape[0]
        new_meta['NAXIS2']=new_data.shape[1]
        new_meta['NAXIS3']=new_data.shape[2]
        return Cube(new_data,new_meta)

    def featformat(self):
        adn=np.meshgrid(self.ra_axis,self.dec_axis,self.nu_axis,indexing='ij')
        X=np.empty(3,len(self_ra_axis)*len(self.dec_axis)*len(self.nu_axis))
        X[0]=adn[0].flatten()
        X[1]=adn[1].flatten()
        X[2]=adn[2].flatten()
        y=self.data.flatten()
        return X,y


#   def add(self,function,center,params,window):
       
