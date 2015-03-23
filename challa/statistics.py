import numpy as np

#def gen_line(spe_form, freq, axis):
#    """ 
#    Returns a spectral line distribution and its application window.
#
#    gen_line() generates a normalized distribution (with compact support) 
#    of a spectral line centered at freq, within an axis. The spe_form 
#    parameter is a tuple where the first element is the name of the 
#    distribution, and the next elements are their parameters. 
#    
#    The available distributions are:
#    - ("skew",fwhm,alpha) : a skew-normal distribution, where fwhm is 
#              the full width at half maximum, and alpha is the curtosis 
#              parameter. If alpha = 0, it degenerates to a Gaussian 
#              distribution, if alpha < 0, the distribution is left-biased
#              and alpha > 0 means a right bias.
#    - TODO: implement more spectral distributions
#    
#    The function returns the tuple (distro,window), where distro is the
#    spectral distribution (array), and window is a tuple of upper and 
#    lower indices (see freq_window function). The window factor is 3*sigma
#    to include relatively large tails.
#    @rtype : tuple
#    @param spe_form: a spatial form (see doc)
#    @param freq: central frequency
#    @param axis: a equispaced array of elements
#    """
#
#    def pdf(x):
#        return 1 / sqrt(2 * pi) * exp(-x ** 2 / 2)
#
#    def cdf(x):
#        return (1 + erf(x / sqrt(2))) / 2
#
#    def skew(x, ep=0, wp=1, ap=0):
#        t = (x - ep) / wp
#        return 2 / wp * pdf(t) * cdf(ap * t)
#
#    fwhm = spe_form[1]
#    a = spe_form[2]
#
#    sigma = fwhm2sigma(freq,fwhm)
#    factor = 3 * sigma
#    window = freq_window(freq, factor, axis)
#    if window[0] > window[1]:
#        return False, window
#
#    d = a / sqrt(1.0 + a ** 2)
#    w = sigma / sqrt(1.0 - (2.0 / pi) * d ** 2)
#    e = freq - w * d * sqrt(2.0 / pi)
#    distro = skew(axis[window[0]:window[1] + 1], e, w, a)
#    ss = sum(distro)
#    if ss != 0:
#        distro = distro / ss
#    #distro=np.sqrt(2*np.pi)*sigma*distro
#    return distro, window


class Gaussian:
    def __init__(self,(a,b,mu,Sigma),inv_sigma=False):
       if inv_sigma:
          self.Sigma=False
          self.Lambda=Sigma
       else:
          self.Sigma=Sigma
          self.Lambda=np.linalg.inv(Sigma)
       self.a=a
       self.b=b
       self.mu=mu


    def getCovariance(self):
       if self.Sigma==False:
          self.Sigma=np.linalg.inv(Lambda)
       return Sigma

    def getPrecision(self):
       return self.Lambda

    def evaluate(self,X,consider_b=True):
       C=np.empty_like(X)
       C[0]=X[0] - self.mu[0]
       C[1]=X[1] - self.mu[1]
       C[2]=X[2] - self.mu[2]
       V=C*(self.Lambda.dot(C))
       quad=V.sum(axis=0)
       #print "quad"
       #print quad
       v=np.exp(-quad/2)
       retval=self.a*v
       #print "retval"
       #print retval
       if consider_b:
          retval= retval + self.b
       return retval

