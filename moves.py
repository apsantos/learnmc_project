import math as ma
import numpy as np
import random as rand
class moves(object):
  def __init__(self,n,box,energy,temp, mu):
    """Initiate moves class

    Parameters
    ----------
    n : int
        number of molecules
    box : class
       simulation box parameter information, see box.py
    energy : class
       simulation energy information, see energy.py
    move : class
       move routine and information, see moves.py
    temp : float
        temperature
    mu : float
        chemical potential

    self.nattempt : 1d array 
        number of attempted moves of each type
    self.naccept : 1d array 
        number of accepted moves of each type

    """
    self.nmol = n
    self.box = box
    self.ndim = box.ndim
    self.vol = box.vol
    self.rho = self.nmol / self.vol
    self.energy = energy
    self.beta = 1./float(temp)
    self.mu = mu
    self.expbetamu = ma.exp(self.beta*self.mu)
    self.nmovetypes = 3
    self.nattempt = np.zeros((self.nmovetypes))
    self.naccept = np.zeros((self.nmovetypes))

  def setmovemix(self,displace,insdel):
    """Set up the move probability and mixes

    Parameters
    ----------
    displace : float
        fraction of MC moves for displacment
    insdel : float
        fraction of MC moves for insertion and deletion (50% split between insertion and deletion)

    """
    self.ifrac_displacement = displace
    self.ifrac_insert = insdel*0.5
    self.ifrac_delete = insdel*0.5
    self.frac_displacement = displace
    self.frac_insert = insdel*0.5 + self.frac_displacement
    self.frac_delete = insdel*0.5 + self.frac_insert
    if (self.frac_delete != 1):
        print('learnmc error: move mix fractions do not add to 1')
        sys.exit(main())    

  def copypos(self, pos):
    """Deep python copy so pos is not changed else where

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues

    """
    npos = np.zeros((len(pos[:,0]),self.ndim))
    for imol in range(self.nmol):
      for idim in range(self.ndim):
        npos[imol,idim] = pos[imol,idim]
    return npos

  def copyexist(self, ex):
    """Deep python copy so the exist array is not changed else where

    Parameters
    ----------
    ex : 1d numpy array (maxnmol length)
        0 if in the resevoir (not in the box) 
        1 if in the box

    """
    nex = np.zeros((len(ex[:])))
    for imol in range(self.nmol):
      nex[imol] = ex[imol]
    return nex

  def move(self, pos, ex):
    """Metropolis algorithm
    0) measure the old energy
    1) pick a random number
    2) move
    3) measure the new energy
    4) apply move-appropriate acceptance criteria
    5) update arrays

    Parameters
    ----------
    pos : 2d numpy array 
        positions of all molecues
    ex : 1d numpy array (maxnmol length)
        0 if in the resevoir (not in the box) 
        1 if in the box

    """
    tfrac = rand.random()
    movetype = 0 # PROJECT

    if movetype != 1: 
      if self.nmol == 0:
        self.energy.tu = 0.0
        self.nattempt[movetype] += 1
        return

      im = rand.randint(0,self.nmol-1)
      oldenergy = self.energy.calcenergy(pos,im,self.nmol)
    else:
      oldenergy = 0.0

    npos = self.copypos(pos)
    nex = self.copyexist(ex)

    if movetype == 0:
      npos[:,:] = self.displace(npos[:,:], im)
      n_nmol = self.nmol
    elif movetype == 1:
      npos[:,:], nex, im = self.insert(npos, nex)
      n_nmol = self.nmol + 1
    elif movetype == 2:
      npos[:,:], nex = self.delete(npos, nex, im)
      n_nmol = self.nmol - 1

    if movetype != 2:
      newenergy = self.energy.calcenergy(npos,im,int(n_nmol))
    else:
      newenergy = 0.0
    # apply metropolis criteria
    try:
      metrop = 0 # PROJECT
    except OverflowError:
      metrop = float('inf')

    if movetype == 1:
      metrop *= 1 # PROJECT

    if movetype == 2:
      metrop *= 1 # PROJECT

    tfrac2 = rand.random()
    if tfrac2 < metrop:
      for idim in range(self.ndim):
        pos[im,idim] = npos[im,idim]
      self.energy.tu = newenergy
      self.naccept[movetype] += 1
      self.nmol = int(n_nmol)
    else:
      self.energy.tu = oldenergy

    self.nattempt[movetype] += 1

  def insert(self,ipos,iexist):
    """Insert a molecule to random location in simulation volume

    Parameters
    ----------
    ipos : 2d numpy array 
        positions of all molecues
    iex : 1d numpy array (maxnmol length)
        0 if in the resevoir (not in the box) 
        1 if in the box

    Returns
    -------
    ipos : 2d numpy array 
        new positions of all molecues
    iex : 1d numpy array (maxnmol length)
        new version
    tnm : int
        which molecule was inserted

    """
    tnm = self.nmol
    for idim in range(self.box.ndim):
      ipos[tnm,idim] = rand.uniform(self.box.box[0,idim],self.box.box[1,idim])
    iexist[tnm] = 1
    return ipos, iexist, tnm

  def delete(self,delpos,delexist,im):
    """Delete a random molecule 
       Move the last molecule in the list to replace the deleted one

    Parameters
    ----------
    delpos : 2d numpy array 
        positions of all molecues
    delexist : 1d numpy array (maxnmol length)
        0 if in the resevoir (not in the box) 
        1 if in the box
    im : int
        which molecule to delete

    Returns
    -------
    delpos : 2d numpy array 
        new positions of all molecues
    delexist : 1d numpy array (maxnmol length)
        new version

    """
    tnm = self.nmol-1
    # move the one from the back to the deleted
    for idim in range(self.box.ndim):
      delpos[im,idim] = delpos[tnm,idim]

    delexist[tnm] = 0
    return delpos, delexist

  def displace(self,dpos,im):
    """Move molecule to random location in simulation volume

    Parameters
    ----------
    dpos : 2d numpy array 
        positions of all molecues
    im : int
        which molecule to move

    Returns
    -------
    dpos : 2d numpy array 
        new positions of all molecues

    """
    for idim in range(self.box.ndim):
      dpos[im,idim] = rand.uniform(self.box.box[0,idim],self.box.box[1,idim])
    return dpos

  def outputefficiency(self, fname='mc'):
    """Output the efficiency of each move

    Parameters
    ----------
    fname : string
        out filename beginning

    """
    movefile = open(fname + '.moves', 'w')
    movefile.write('# fa_translate fa_insert fa_delete\n')
    for imovetype in range(self.nmovetypes):    
      if self.nattempt[imovetype] > 0:
        ratio = self.naccept[imovetype]/float(self.nattempt[imovetype])
      else:
        ratio = 0
      movefile.write('%f ' % ratio)
    movefile.write('\n')
