import random
import numpy as np
import pylab as plt
import matplotlib.animation as animation
from astropy.io import fits
from astropy  import constants as const
from scipy import pi, sqrt, exp
import math
import db
from scipy.special import erf
from astropy import units as u
from common import *

# ## ALMA Specific Constants ###
MAX_CHANNELS = 9000
MAX_BW = 2000000000
ALMA_bands = {'3': [88000, 116000], '4': [125000, 163000], '6': [211000, 275000], '7': [275000, 373000],
              '8': [385000, 500000], '9': [602000, 720000]}
ALMA_noises = {'3': 0.01, '4': 0.012, '6': 0.02, '7': 0.04, '8': 0.08, '9': 0.16}

### IMC Specific Constants ###
inten_group = [('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')]
inten_values = [[0.1, 2], [20, 60], [5, 20], [1, 10]]
default_iso_abundance = {'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': 1.0 / 120, '34S': 1.0 / 30, '33S': 1.0 / 120,
                         '13N': 1.0 / 30, 'D': 1.0 / 30}

### Core ASYDO Classes ###

class Universe:
    """
    A synthetic universe where to put synthetic objects.
    """

    def __init__(self, log):
        """
        The log parameter is an opened file descriptor for logging purposes
        """
        self.log = log
        self.sources = dict()
    
    @u.quantity_input(ra=u.deg)
    @u.quantity_input(dec=u.deg)
    def create_source(self, name, ra, dec):
        """
        A source needs a name and a spatial position (ra,dec).
        """
        self.sources[name] = _Source(self.log, name, ra, dec)

    def add_component(self, source_name, model):
        """
        To add a component a Component object must be instantiated (model), and added to 
        a source called source_name.
        """
        self.sources[source_name].add_component(model)

    @u.quantity_input(pos_ra=u.deg)
    @u.quantity_input(pos_dec=u.deg)
    @u.quantity_input(pos_nu=u.MHz)
    @u.quantity_input(delta_ra=u.deg/u.pix)
    @u.quantity_input(delta_dec=u.deg/u.pix)
    @u.quantity_input(delta_nu=u.MHz/u.pix)
    @u.quantity_input(pix_ra=u.pix)
    @u.quantity_input(pix_dec=u.pix)
    @u.quantity_input(pix_nu=u.pix)
    @u.quantity_input(crpix_ra=u.pix)
    @u.quantity_input(crpix_dec=u.pix)
    @u.quantity_input(crpix_nu=u.pix)
    def gen_cube(self, name, pos_ra,pos_dec,pos_nu,delta_ra,delta_dec,delta_nu,pix_ra,pix_dec,pix_nu,crpix_ra,crpix_dec,crpix_nu):
        """
        Returns a SpectralCube object where all the sources within the FOV and BW are projected.

        This function needs the following parameters:
        - name    : name of the cube
        - pos     : position (nu,dec,ra)
        - pix     : pixels (nu,dec,ra)
        - res     : resolution (nu,dec,ra)
        - crpix   : center pixed (nu,dec,ra) from 1
        """
        cube = _synthetic_cube(self.log, name, pos,res,pix,crpix)
        for src in self.sources:
            self.log.write('*** Source: ' + src + '\n')
            self.sources[src].project(cube)
        return cube

    def save_cube(self, cube, filename):
        """
        Wrapper function that saves a cube into a FITS (filename).
        """
        self.log.write('   -++ Saving FITS: ' + filename + '\n')
        cube.save_fits(filename)

    def remove_source(self, name):
        """
        Deletes a source and its components.
        """
        self.log.write('Removing source ' + name)
        return self.sources.remove(name)


class _Source:
    """
    A generic source of electromagnetic waves with several components.
    """

    def __init__(self, log, name, alpha, delta):
        """ Parameters:
               * log: logging descriptor
               * name: a name of the source
               * alpha: right ascension 
               * delta: declination
        """
        self.log = log
        log.write('+++ Source \'' + name + '\' added\n')
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def add_component(self, model):
        """ Defines a new component from a model.
        """
        code = self.name + '-c' + str(len(self.comp) + 1)  #+ '-r' + str(self.alpha) +'-d'+str(self.delta)
        self.comp.append(model)
        model.register(code, self.alpha, self.delta)
        self.log.write(' |- Component added to \'' + self.name + '\' (' + code + ')\n')
        self.log.write(' ---+ Model: ' + model.info() + '\n')

    def project(self, cube):
        """
        Projects all components in the source to a cube.
        """
        for component in self.comp:
            self.log.write('  |- Projecting ' + component.comp_name + '\n')
            component.project(cube);



def _synthetic_cube(log, name,pos,res,pix,crpix, band_freq=ALMA_bands,
                 band_noises=ALMA_noises):
        """ 
        Obligatory Parameters:
        - log	  : descriptor of a log file
        - name    : name of the cube 
        - pos     : position (ra,dec,nu)
        - pix     : pixels (ra,dec,nu)
        - res     : resolution (ra,dec,nu)
        - crpix   : center pixed (ra,dec,nu) from 1
                
        Optional Parameters:
        - band_freq   : a diccionary of frequency ranges for the bands (key = band_name, value = (lower,upper)) (IN MHZ!!!)
        - band_noises : a dictionary of noise levels for each band (key = band_name, value = noise)
        """
        log.write('[*] Generating Synthetic Cube ' + name + '\n')
        log.write('  |- Coordinates : ra=' + str(pos[0]) + '[deg] dec=' + str(pos[1]) + '[deg] nu='+str(pos[2])+'[Hz] \n')
        if alpha > 90 or alpha < -90:
            raise Exception('!!! ERROR: invalid coordinate: ra=' + alpha)
        if delta > 90 or delta < -90:
            raise Exception('!!! ERROR: invalid coordinate: dec=' + delta)
        [fu,fd]=[pos[2]+(pix[2]-(crpix[2]-1))*res[2],pos[2]- (crpix[2]-1)*res[2]]
        [du,dd]=[pos[1]+(pix[1]-(crpix[1]-1))*res[1],pos[1]- (crpix[1]-1)*res[1]]
        [au,ad]=[pos[0]+(pix[2]-(crpix[0]-1))*res[0],pos[0]- (crpix[0]-1)*res[0]]
        bw=fu - fd;
        log.write('  |- Bandwidth : bw=' + str(bw))
        if  bw  > MAX_BW:
            log.write('!!! WARNING: max ALMA bandwidth exceeded\n')
        if  pix[0] > MAX_CHANNELS:
            log.write('!!! WARNING: max ALMA channels exceeded\n')
        band = 'NO_BAND'
        for bnd in band_freq:
            freqs = band_freq[bnd]
            if fd >= freqs[0] and fu <= freqs[1]:
                band = bnd
                log.write('  |- Band: ' + bnd + '\n')
        if band == 'NO_BAND':
            log.write('!!! WARNING: not in a valid ALMA band\n')
        if band == 'NO_BAND':
            noise = 0.0001
        else:
            noise = band_noises[self.band]
        data = (np.random.random(pix) - 0.5 * np.ones(pix))* 2 * noise
        prihdr = fits.Header()
        prihdr['AUTHOR'] = 'Astronomical SYnthetic Data Observatory'
        prihdr['COMMENT'] = name
        prihdr['SIMPLE'] = True
        # prihdr['BITPIX'] = 8
        prihdr['NAXIS'] = 4
        prihdr['NAXIS1'] = pix[2]
        prihdr['NAXIS2'] = pix[1]
        prihdr['NAXIS3'] = pix[0]
        prihdr['BMAJ'] = res[1]
        prihdr['BMIN'] = res[1]
        prihdr['CTYPE1'] = 'RA---SIN'
        prihdr['CRVAL1'] = pos[0]
        prihdr['CDELT1'] = res[0]
        prihdr['CRPIX1'] = crpix[0]
        prihdr['CUNIT1'] = 'deg'
        prihdr['CTYPE2'] = 'DEC--SIN'
        prihdr['CRVAL2'] = pos[1]
        prihdr['CDELT2'] = res[1]
        prihdr['CRPIX2'] = crpix[1]
        prihdr['CUNIT2'] = 'deg'
        prihdr['CTYPE3'] = 'FREQ'
        prihdr['CRVAL3'] = pos[2] 
        prihdr['CDELT3'] = res[2]
        prihdr['CRPIX3'] = crpix[2]
        prihdr['CUNIT3'] = 'Hz'
        prihdr['CTYPE4'] = 'STOKES'
        prihdr['CRVAL4'] = 1.0
        prihdr['CDELT4'] = 1.0
        prihdr['CRPIX4'] = 1.0
        return spaectral.Cube(data,prihdr)

class Component:
    """Abstract component model"""

    def __init__(self, log, z_base=0.0):
        """ log: file descriptor for logging
            z_base[=0] : optional parameter to set base redshift (if not, please use set_radial_vel)
        """
        self.log = log
        self.z = z_base
        self.rv = const.c.value * ((self.z ** 2 + 2 * self.z) / (self.z ** 2 + 2 * self.z + 2))

    def set_radial_velocity(self, rvel):
        """Set radial velocity rvel in m/s"""
        self.rv = rvel
        self.z = math.sqrt((1 + self.rv/const.c.value) / (1 - self.rv / const.c.value)) - 1

    def info(self):
        """Print relevant information of the component"""
        return "(none)"

    def register(self, comp_name, alpha, delta):
        """Register the component name and angular position (alpha,delta)"""
        self.comp_name = comp_name
        self.alpha = alpha
        self.delta = delta

    def project(self, cube):
        """Project the component in the cube"""
        pass


class IMCM(Component):
    """ Interstellar Molecular Cloud Model """

    def __init__(self, log, dbpath, mol_list, temp, form, z_grad, z_base=0.0, abun_max=10 ** -5,
                 abun_min=10 ** -6, abun_CO=1.0, iso_abun=default_iso_abundance):
        Component.__init__(self, log, z_base)
        self.form = form
        self.z_grad = z_grad
        self.dbpath = dbpath
        self.temp = temp
        self.intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(abun_min, abun_max)
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
                abun += abun_CO
            for iso in iso_abun:
                if iso in mol:
                    abun *= iso_abun[iso]
            self.intens[mol] = abun

    def change_intensities(self, intens):
        '''User defined dictionary in the form {molecule: intensity}'''
        self.intens = intens;

    def info(self):
        return "mol_list = " + str(self.intens.keys()) + " @ form=" + str(self.form) + ", z=" + str(self.z) + ", grad=" + str(self.z_grad)

    def project(self, cube):
        arr_code = []
        arr_mol = []
        arr_chname = []
        arr_rest_freq = []
        arr_rad_vel = []
        arr_fwhm = []
        arr_temp = []
        #self.log.write('   --+ Generating template image\n')  # TODO More info
        #T, (ybord, xbord) = gen_surface(self.spa_form, self.alpha, self.delta, cube.alpha_axis, cube.delta_axis)
        #if isinstance(T, bool):
        #    return
        #G = gen_gradient(self.z_grad, self.alpha, self.delta, cube.alpha_axis, cube.delta_axis, ybord, xbord)
        #self.log.write('   --+ Generating template line distribution\n')  #TODO More info
        #self.log.write('   --+ Loading and correcting lines\n')
        dba = db.lineDB(self.dbpath)
        dba.connect()
        freq_init_corr = cube.nu_axis[0] / (1 + self.z)
        freq_end_corr = cube.nu_axis[-1] / (1 + self.z)
        counter = 0
        used = False
        for mol in self.intens:
            # For each molecule specified in the dictionary
            # load its spectral lines
            linlist = dba.getSpeciesLines(mol, freq_init_corr,freq_end_corr)  # Selected spectral lines for this molecule
            rinte = inten_values[0]
            for j in range(len(inten_group)):  # TODO baaad python...
                if mol in inten_group[j]:
                    rinte = inten_values[j]
            rinte = random.uniform(rinte[0], rinte[1])

            for lin in linlist:
                counter += 1
                trans_temp = lin[5]
                temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                if temp < cube.noise:
                    continue
                freq = (1 + self.z) * lin[3]*1000000.0  # Catalogs are in Mhz
                self.log.write('      |- Projecting ' + str(lin[2]) + ' (' + str(lin[1]) + ') around ' + str(
                    freq/1000000) + ' Mhz, at ' + str(temp) + ' K\n')
                (X,(n0,n1,d0,d1,r0,r1)) = cube.feature_space(xmax,2*params['weight_deltas']*resv)
                par=[temp,0,self.alpha,self.delta,freq,self.form[3],form[0],form[1],form[2],z_grad[0],z_grad[1]]
                G=Gaussian(clumps.to_gauss(par),True)
                M=G.evaluate(X,False).reshape((n1-n0+1,d1-d0+1,r1-r0+1))
                cube.data[n0:n1+1,d0:d1+1,r0:r1+1] += M
                #for xp in range(xbord[0], xbord[1]):
                #    for yp in range(ybord[0], ybord[1]):
                #        freq = spectral.doppler(lin[3],self.rv + G[yp - ybord[0], xp - xbord[0]])
                #        L, Lbord = gen_line(self.spe_form, freq, cube.freq_axis)
                #        if isinstance(L, bool):
                #            continue
                #        cube.data[Lbord[0]:Lbord[1] + 1, yp, xp] = cube.data[Lbord[0]:Lbord[1] + 1, yp, xp] + T[yp - ybord[0], xp - xbord[0]] * temp* L
                used = True
                arr_code.append(self.comp_name + '-r' + str(self.alpha) + '-d' + str(self.delta) + "-l" + str(counter))
                arr_mol.append(mol)
                arr_temp.append(temp)
                arr_chname.append(str(lin[2]))
                arr_rest_freq.append(str(lin[3]))
                arr_rad_vel.append(self.rv)
                arr_fwhm.append(self.spe_form[1])
        dba.disconnect()
        if not used:
            return
        hduT = fits.PrimaryHDU()
        hduT.data = T;
        hduG = fits.PrimaryHDU()
        hduG.data = G;
        tbhdu = fits.new_table(fits.ColDefs([
            fits.Column(name='line_code', format='60A', array=arr_code),
            fits.Column(name='mol', format='20A', array=arr_mol), \
            fits.Column(name='chname', format='40A', array=arr_chname), \
            fits.Column(name='rest_freq', format='D', array=arr_rest_freq), \
            fits.Column(name='rad_vel', format='D', array=arr_rad_vel), \
            fits.Column(name='fwhm', format='D', array=arr_fwhm), \
            fits.Column(name='temp', format='D', array=arr_temp)]))
        cube._add_HDU(hduT)
        cube._add_HDU(hduG)
        cube._add_HDU(tbhdu)



