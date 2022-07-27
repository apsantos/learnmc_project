import numpy as np
import random as rand
import sys
import matplotlib.pyplot as plt

class movie(object):
  def __init__(self,s,box):
    """Initiate Movie class

    Parameters
    ----------
    s : float
        the lennard-jones sigma
    box : class
       simulation box parameter information, see box.py

    """
    yfactor = (box.box[1,1]-box.box[0,1])/(box.box[1,0]-box.box[0,0])
    self.size = (s*23*10/(box.box[1,0]-box.box[0,0]))**2.
    self.box = box.box

    self.figsize = (5,5*yfactor)

  def update(self,pos,exist):
    """Update the movie

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    exist : 1d numpy array
       either 1 or 0 (not existant), used for muVT

    """
    # ims is a list of lists, each row is a list of artists to draw in the
    # current frame; here we are just animating one artist, the image, in
    # each frame
    plt.subplots(figsize=self.figsize)
    tnmol = sum(exist)
    im = plt.scatter(pos[:tnmol,0], pos[:tnmol,1], s=self.size)
    plt.xlim(self.box[:,0])
    plt.ylim(self.box[:,1])
    plt.show()
        
class initialize(object):
  def __init__(self,maxn):
    self.maxn = maxn

  def random(self,nmol,box):
    """Generate random initial configuration

    Parameters
    ----------
    nmol : int
        number of molecules
    box : class
       simulation box parameter information, see box.py

    """
    if nmol > self.maxn:
        print('ERROR: Max number of molecules is: %d' % self.maxn)
        sys.exit()

    pos = np.zeros((self.maxn,box.ndim))
    posexist = np.zeros((self.maxn))
    for imol in range(nmol):
      for idim in range(box.ndim):
        pos[imol,idim] = rand.uniform(box.box[0,idim],box.box[1,idim])
        posexist[imol] = 1
    return pos, posexist

class configfile(object):
  def __init__(self,ndim, fname='config'):
    """Initiate configuration file class

    Parameters
    ----------
    ndim : int
        number of simulation dimensions
    fname : string
       file output name

    """
    self.file = open(fname + '.xyz', 'w')
    self.ndim = ndim

  def writeframe(self,imove,pos, nmol):
    """Generate random initial configuration

    Parameters
    ----------
    imove : int
        MC move number
    pos : 2d numpy array 
        positions of all molecues
    nmol : int
        number of molecules

    """
    self.file.write('%d\nframe %d\n' % (nmol, imove))
    if self.ndim == 3:
      for imol in range(nmol):
        self.file.write('C %f %f %f\n' % (pos[imol,0],pos[imol,1],pos[imol,2]))
    elif self.ndim == 2:
      for imol in range(nmol):
        self.file.write('C %f %f 0.0\n' % (pos[imol,0],pos[imol,1]))
    elif self.ndim == 1:
      for imol in range(nmol):
        self.file.write('C %f 0.0 0.0\n' % (pos[imol,0]))
