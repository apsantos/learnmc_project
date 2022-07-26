import math as ma
import numpy as np
class energy(object):
  def __init__(self,nmol,box, temp):
    """Initiate energy class
    Also used for pressure and virial

    Parameters
    ----------
    nmol : int
        number of molecules
    box : class
       simulation box parameter information, see box.py
    temp : float
        temperature

    """
    yfactor = (box.box[1,1]-box.box[0,1])/(box.box[1,0]-box.box[0,0])
    self.calctensor = False
    self.nmol = nmol
    self.box = box
    self.ndim = len(box.box[0,:])
    self.vol = box.vol
    self.rho = nmol/box.vol
    self.T = temp

  def setparameters(self,sigma=1.,epsilon=1.,cutoff=-1):    
    """Set up the interaction parameters

    Parameters
    ----------
    sigma : float
        LJ distance
    epsilon : float
        LJ interaction strength
    cutoff : float
        interaction cutoff distance

    """
    self.a = 4.*epsilon*sigma**12.
    self.b = 4.*epsilon*sigma**6.
    self.fprea = 48. * epsilon * sigma**12.
    self.fpreb = 24. * epsilon * sigma**6.
    if cutoff < 0:
      self.rc =(self.box[1,0]/2.)
    else:
      self.rc = cutoff
    self.rc2 = self.rc * self.rc
    # if you want to add tail interactions
    # self.rho = self.nmol/self.vol
    # factor = 8./3.*3.14159*self.rho
    # etail = factor*((1/3.)*self.rc**9.-self.rc**3.)
    # ptail = 2.*factor *self.rho*((2./3.)*self.rc**9.-self.rc**3.)

  def calcforces(self,pos,nmol):
    """Calculate the forces

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    nmol : int
        number of molecules

    """
    self.f = np.zeros((nmol,self.ndim))
    r = np.zeros((3))
    for imol in range(0,nmol-1):
      for jmol in range(imol+1,nmol):
        for idim in range(self.ndim):
          r[idim] = pos[imol,idim] - pos[jmol,idim]
        self.box.applypbc(r)
        r2 = 0.0
        for idim in range(self.ndim):
          r2 += r[idim]**2.0
        if r2 < self.rc2:
          r6 = r2**3.
          r12 = r6*r6
          tf = (self.fprea/r12 - self.fpreb/r6)/r2*r
          self.f[imol,:] += tf
          self.f[jmol,:] -= tf
          self.vir += np.dot(r, tf)
          if self.calctensor:
            for idim in range(self.ndim):
              for jdim in range(self.ndim):
                self.virtensor[idim,jdim] += r[idim] * tf[jdim]

  def calcvirial(self, pos, nmol):
    """Calculate the virial pressure

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    nmol : int
        number of molecules

    """
    self.vir = 0.0
    self.virtensor = np.zeros((self.ndim, self.ndim))
    self.calcforces(pos,nmol)
    self.vir /= float(self.ndim)

  def calcpressure(self, pos, nmol, vol):
    """Calculate the pressure

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    nmol : int
        number of molecules
    vol : float
        simulation box volume

    """
    self.calcvirial(pos, nmol)
    rho = nmol/vol
    p = self.T * rho + self.vir/vol
    self.tp = p

  def calcfullenergy(self,pos,nmol):
    """Calculate the system, looping over all molecules

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    nmol : int
        number of molecules

    """
    self.tu = 0
    r = np.zeros((3))
    for imol in range(0,nmol-1):
      for jmol in range(imol+1,nmol):
        for idim in range(self.ndim):
          r[idim] = pos[imol,idim] - pos[jmol,idim]
        self.box.applypbc(r)
        r2 = 0.0
        for idim in range(self.ndim):
          r2 += r[idim]**2.0

        if r2 < self.rc2:
          r6 = r2**3.
          r12 = r6*r6
          self.tu += self.a/r12 - self.b/r6

  def calcenergy(self,pos,mmol,nmol):
    """Calculate the energy of test molecules

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    mmol : int
        test molecule number
    nmol : int
        number of molecules

    """
    mollist = [*range(0,mmol),*range(mmol+1,nmol)]
    mpos = pos[mmol,:]
    e = 0.0
    r = np.zeros((3))
    for imol in mollist:
      for idim in range(self.ndim):
        r[idim] = pos[imol,idim] - mpos[idim]
      self.box.applypbc(r)
      r2 = 0.0
      for idim in range(self.ndim):
        r2 += r[idim]**2.0

      if r2 < self.rc2:
        r6 = r2**3.
        r12 = r6*r6
        e += self.a/r12 - self.b/r6

    return e
