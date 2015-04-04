import numpy as np
import copy

from astropy import constants as const
from astropy import units as u
from astropy.io import fits 

# ## Helper constants ###
#SPEED_OF_LIGHT = 299792458.0
S_FACTOR = 2.354820045031  # sqrt(8*ln2)
#DEG2ARCSEC = 3600.0

def fwhm2sigma(freq,fwhm):
    """
    Compute the sigma in Hz given a frequency in Hz and a fwhm in m/s
    """
    sigma = (fwhm / S_FACTOR) * (freq / const.c.value)
    return sigma

def doppler(freq,rv):
    freq_new = math.sqrt((1 + rv / const.c.value) / (1 - rv / const.c.value)) * freq
    return freq_new



#def gen_gradient(form, alpha, delta, alpha_axis, delta_axis, ybord, xbord):
#    """
#    Returns a matrix with gradient values within the axes, restricted to y_bord and x_bord.
#
#    gen_gradient() generate a gradient matrix centered at (alpha,delta), for the alpha_axis
#    and delta_axis, restricted by the windows ybord and xbord. The form parameter is a tuple 
#    where the first element is the gradient function, and the next elements are their parameters. 
#    The available gradient functions are:
#    - ("linear", theta, m) : linear gradient, rotated by theta, and intensity m (km/s/arcsec)
#    """
#    np.set_printoptions(threshold=np.nan)
#    gtype = form[0]
#    theta = form[1]
#    km_sarcs = form[2]
#    alpha_axis = alpha_axis[xbord[0]:xbord[1]]
#    delta_axis = delta_axis[ybord[0]:ybord[1]]
#    alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse=False, indexing='xy')
#    Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
#    Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
#    XX = Xc * math.cos(-theta) - (Yc) * math.sin(-theta);
#    res = np.reshape(km_sarcs * XX, (len(delta_axis), len(alpha_axis)))
#    return res



#def gen_surface(form, alpha, delta, alpha_axis, delta_axis):
#    """ 
#    Returns a spatial distribution and its application window.
#
#    gen_surface() generates a standarized surface distribution 
#    (with compact support) centered at (alpha,delta), within the axes. 
#    The spa_form parameter is a tuple where the first element is the 
#    name of the distribution, and the next elements are their parameters. 
#    
#    The available distributions are:
#    - ("normal",sx,sy,theta) : a 2D Gaussian (normal) distribution, 
#                with diagonal covariance [[sx,0],[0,sy], rotated by 
#                theta (in radians!).
#    - ("exp",sx,sy,theta) : exponential distribution, with semiaxes
#                sx and sy, rotated by theta (in radians!).
#    - TODO: implement more spatial distributions
#                               
#    The function returns the tuple (distro,ybord,xbord), where distro is 
#    the spatial distribution (matrix), and ybord and xbord are the tuples 
#    of upper and lower indices for the axes (see spatial_window function).
#    If the window is spurious (out of the axes or too small) this function
#    returns False.   
#    """
#    stype = form[0]
#    sx = form[1] / DEG2ARCSEC
#    sy = form[2] / DEG2ARCSEC
#    theta = form[3]
#    spx = abs(3 * sx * math.cos(theta)) + abs(3 * sy * math.sin(theta))
#    spy = abs(3 * sx * math.sin(theta)) + abs(3 * sy * math.cos(theta))
#    ybord, xbord = spatial_window(alpha, delta, spx, spy, alpha_axis, delta_axis)
#    if xbord[0] > xbord[1] or ybord[0] > ybord[1]:
#        return False, [ybord, xbord]
#    alpha_axis = alpha_axis[xbord[0]:xbord[1] + 1]
#    delta_axis = delta_axis[ybord[0]:ybord[1] + 1]
#    alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse=False, indexing='xy')
#    Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
#    Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
#    XX = (Xc) * math.cos(-theta) - (Yc) * math.sin(-theta)
#    YY = (Xc) * math.sin(-theta) + (Yc) * math.cos(-theta)
#        u = (XX / sx) ** 2 + (YY / sy) ** 2
#        sol = sx * sy * np.exp(-u / 2) / (2 * math.pi)
#    elif stype == 'exp':
#        u = sqrt((XX / sx) ** 2 + (YY / sy) ** 2)
#        sol = sx * sy * np.exp(-u / sqrt(2))
#    else:
#        print('!!! ERROR: No such surface type')
#        return False, [ybord, xbord]
#    mm = max(sol)
#    if mm != 0:
#        sol = sol / mm
#    res = np.reshape(sol, (len(delta_axis), len(alpha_axis)))
#    return res, (ybord, xbord)




#def spatial_window(alpha, delta, spx, spy, alpha_axis, delta_axis):
#    """
#    Compute a 2D window centered at (alpha,delta) within alpha_axis and delta_axis.
#
#    spatial_window() returns a tuple (alpha_window,delta_window), where each element 
#    are the lower and upper indices for each axis. The window correspond to a window
#    of freq +- spx or freq +- spy depending on the axis. This window is limited by 
#    the axis borders. Equispaced axes are assumed.
#    @rtype : tuple
#    @param alpha: RA center
#    @param delta: DEC center
#    @param spx: Semiaxis X
#    @param spy: Semiaxis Y
#    @param alpha_axis: a equispaced array of elements
#    @param delta_axis: a equispaced array of elements
#    """
#    dlta_x = alpha_axis[1] - alpha_axis[0]
#    dlta_y = delta_axis[1] - delta_axis[0]
#    xbord = [int(round((alpha - spx - alpha_axis[0]) / dlta_x)), int(round((alpha + spx - alpha_axis[0]) / dlta_x))]
#    ybord = [int(round((delta - spy - delta_axis[0]) / dlta_y)), int(round((delta + spy - delta_axis[0]) / dlta_y))]
#    if xbord[0] < 0:
#        xbord[0] = 0
#    if ybord[0] < 0:
#        ybord[0] = 0
#    if xbord[1] > len(alpha_axis):
#        xbord[1] = len(alpha_axis)
#    if ybord[1] > len(delta_axis):
#        ybord[1] = len(delta_axis)
#    return ybord, xbord


#def freq_window(freq, factor, axis):
#    """
#    Compute a window centered at freq within an axis.
#
#    freq_window() returns a tuple with the lower and upper indices in the axis for a
#    frequency window of freq +- factor. The size is at most 2*factor, but
#    is limited by the axis borders. It returns (0,0) or (end,end) if the window is
#    out of the axis by the left or right respectively (i.e., end = len(axis)).
#    Equispaced axes are assumed.
#    @rtype : tuple
#    @param freq: center frequency
#    @param factor: the window factor (half window size)
#    @param axis: a equispaced array of elements
#    """
#    dlta = axis[1] - axis[0]
#    ini = int(round((freq - factor - axis[0]) / dlta))
#    end = int(round((freq + factor - axis[0]) / dlta))
#    if ini < 0:
#        ini = 0
#    if end > len(axis):
#        end = len(axis)
#    return ini, end

def cube_data_unravel(data,idx):
   return data.reshape((idx[5]-idx[4]+1,idx[3]-idx[2]+1,idx[1]-idx[0]+1))

def cube_data_stack(data):
   return data.sum(axis=0)


class Cube:
    """
    A spectral cube.
    """
    def __init__(self,data,meta):
        """ data = numpy 3d array 
            meta = header of fits
        """
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
        hdu = fits.PrimaryHDU(header=self.meta)
        hdu.data = self.data
        self.hdulist = fits.HDUList([hdu])
    
    def ravel(self,idx=np.array([])):
        if len(idx)!=6:
           lss=self.data
        else:
           lss=self.data[idx[4]:idx[5]+1,idx[2]:idx[3]+1,idx[0]:idx[1]+1]
        return lss.ravel()

    def stack(self,idx=np.array([])):
        if len(idx)!=6:
           lss=self.data
        else:
           lss=self.data[idx[4]:idx[5]+1,idx[2]:idx[3]+1,idx[0]:idx[1]+1]
        return lss.sum(axis=0)

    def add(self,sc,idx=np.array([])):
        if (len(idx)!=6):
           self.data=self.data + sc
        else:
           self.data[idx[4]:idx[5]+1,idx[2]:idx[3]+1,idx[0]:idx[1]+1] =self.data[idx[4]:idx[5]+1,idx[2]:idx[3]+1,idx[0]:idx[1]+1] + sc
        

    def max(self):
        index=np.unravel_index(self.data.argmax(),self.data.shape)
        y=self.data[index]
        x=np.empty(3)
        x[2]=self.nu_axis[index[0]]
        x[1]=self.dec_axis[index[1]]
        x[0]=self.ra_axis[index[2]]
        return (y,x)
 
    def feature_space(self,center,window):
        ra_ci=np.argmin(np.abs(self.ra_axis-center[0]));
        ra_ui=np.argmin(np.abs(self.ra_axis-center[0]-window[0]));
        ra_li=np.argmin(np.abs(self.ra_axis-center[0]+window[0]));
        dec_ci=np.argmin(np.abs(self.dec_axis-center[1]));
        dec_ui=np.argmin(np.abs(self.dec_axis-center[1]-window[1]));
        dec_li=np.argmin(np.abs(self.dec_axis-center[1]+window[1]));
        nu_ci=np.argmin(np.abs(self.nu_axis-center[2]));
        nu_ui=np.argmin(np.abs(self.nu_axis-center[2]-window[2]));
        nu_li=np.argmin(np.abs(self.nu_axis-center[2]+window[2]));
        

        crval1=self.ra_axis[ra_ci]
        crval2=self.dec_axis[dec_ci]
        crval3=self.nu_axis[nu_ci]
        crpix1=ra_ci - ra_li 
        crpix2=dec_ci - dec_li 
        crpix3=nu_ci - nu_li 
        naxis1=ra_ui-ra_li +1 
        naxis2=dec_ui-dec_li +1
        naxis3=nu_ui-nu_li +1
        ra_axis=np.linspace(crval1-crpix1*self.ra_delta,crval1+(naxis1-crpix1)*self.ra_delta, num=naxis1)
        dec_axis=np.linspace(crval2-crpix2*self.dec_delta,crval2+(naxis2-crpix2)*self.dec_delta, num=naxis2)
        nu_axis=np.linspace (crval3-crpix3*self.nu_delta,crval3+(naxis3-crpix3)*self.nu_delta, num=naxis3)
        adn=np.meshgrid(nu_axis,dec_axis,ra_axis, indexing='ij')
        X=np.empty((3,len(ra_axis)*len(dec_axis)*len(nu_axis)))
        X[2]=adn[0].ravel()
        X[1]=adn[1].ravel()
        X[0]=adn[2].ravel()
        yidx=(ra_li,ra_ui,dec_li,dec_ui,nu_li,nu_ui)
        return X,yidx

    def _add_HDU(self, hdu):
        self.hdulist.append(hdu)

    def save_fits(self, filename):
        """ Simple as that... saves the whole cube """
        self.hdulist.writeto(filename, clobber=True)

    def _updatefig(self, j):
        """ Animate helper function """
        self.im.set_array(self.data[j, :, :])
        return self.im,

    def animate(self, inte, rep=True):
        """ Simple animation of the cube.
            - inte       : time interval between frames
            - rep[=True] : boolean to repeat the animation
          """
        fig = plt.figure()
        self.im = plt.imshow(self.data[0, :, :], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
                             extent=(
                                 self.alpha_border[0], self.alpha_border[1], self.delta_border[0],
                                 self.delta_border[1]))
        ani = animation.FuncAnimation(fig, self._updatefig, frames=range(len(self.freq_axis)), interval=inte, blit=True,
                                      repeat=rep)
        plt.show()

#    def get_spectrum(self, x, y):
#        """ Returns the spectrum of a (x,y) position """
#        xi = int(round((x - self.alpha_axis[0]) / (self.ang_res / DEG2ARCSEC)))
#        yi = int(round((y - self.delta_axis[0]) / (self.ang_res / DEG2ARCSEC)))
#        return self.data[:, yi, xi]

         
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
