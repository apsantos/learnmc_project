import numpy as np
import math as ma
import matplotlib.pyplot as plt
class tally(object):
  def __init__(self,filename='mc'):
    """Initiate tally class
    Used for keeping track of averages also

    Parameters
    ----------
    filename : string
       name of output file

    """
    self.file = open(filename + '.data', 'w')
    self.file.write("# step P U V N\n")
    self.step = []
    self.p = []
    self.u = []
    self.v = []
    self.n = []
    self.rho = []
    self.nsamp = 0

  def average(self):
    """Average and get standard deviation for all properties

    Parameters
    ----------
    self.p : 1d numpy array
       pressure at every sampling point
    self.pave : float
       average pressure thus far
    self.pstd : float
       average pressure standard devation thus far

    """
    self.p = np.array(self.p)
    self.u = np.array(self.u)
    self.v = np.array(self.v)
    self.n = np.array(self.n)
    self.rho = np.array(self.rho)

    self.pave = np.mean(self.p)
    self.uave = np.mean(self.u)
    self.vave = np.mean(self.v)
    self.nave = np.mean(self.n)
    self.rhoave = np.mean(self.rho)

    self.pstd = np.std(self.p)
    self.ustd = np.std(self.u)
    self.vstd = np.std(self.v)
    self.nstd = np.std(self.n)
    self.rhostd = np.std(self.rho)

  def update(self, istep, energy, move):
    """Update the property arrays

    Parameters
    ----------
    istep : int
       this step number
    energy : class
       simulation energy information, see energy.py
    move : class
       move routine and information, see moves.py

    """
    self.nsamp += 1
    P = energy.tp
    U = energy.tu
    V = move.vol
    N = move.nmol
    rho = N / V
    self.step.append(istep)
    self.p.append(P)
    self.u.append(U)
    self.v.append(V)
    self.n.append(N)
    self.rho.append(rho)
    self.file.write('%d %f %f %f %f %f\n' % (istep, P, U, rho, V, N))
    
  def pltseries(self):
    """Plot functions
    """
    plt.plot(self.step, self.p)
    plt.xlabel('step')
    plt.ylabel('P')
    plt.show()

  def outputave(self, T, fname='mc', mu='off'):
    """Write the averages with the standard devations

    Parameters
    ----------
    T : float
       temperature 
    fname : string
        out filename beginning
    mu : either a string or a float
        if we are running a muVT sim this will be a float and the output changes accordingly

    """
    ofile = open(fname + '.ave', 'w')
    ofile.write("# ")
    if mu != 'off':
      ofile.write("mu ")
    ofile.write("T rho P U N V mu_ideal, stds\n")
    if mu != 'off': 
      ofile.write('%f ' % mu)
    ofile.write('%f %f %f %f %f %f %f %f %f %f %f %f %f\n' % (T, self.rhoave, self.pave, self.uave, self.nave, self.vave, -T*ma.log(1./self.rhoave),
                                                           self.rhostd, self.pstd, self.ustd, self.nstd, self.vstd, T*self.rhostd/self.rhoave))
