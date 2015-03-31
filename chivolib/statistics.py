import numpy as np

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
       v=np.exp(-quad/2)
       retval=self.a*v #/v.max()
       if consider_b:
          retval= retval + self.b
       return retval

