import numpy as np
class box(object):
  def __init__(self):
    return
  def setedges(self,edges):
    """Setup the simulation box

    Inputs
    ------
    edges : 1-d numpy array 
       simulation box edges

    Parameters
    ----------
    self.box : 2-d numpy array dimensions (2,ndim)
        the box dimensions self.box[0,0] is the position of the x lower limit
                           self.box[0,1] is the position of the y lower limit
                           self.box[1,0] is the position of the y upper limit
                           self.box[1,1] is the position of the y upper limit
    self.half_boxlength : 1-d numpy array
        half the length in each direction 
    self.boxlength : 1-d numpy array
        length in each direction 
    self.vol : float
        box volume

    """

    self.ndim = len(edges)
    bottomedge = [0] * self.ndim
    self.box = np.array([bottomedge,edges])
    self.boxlength = np.array(edges)
    self.half_boxlength = self.boxlength/2.0
    self.vol = 1.0
    for idim in range(self.ndim):
      self.vol *= (self.box[1,idim]-self.box[0,idim])

  def applypbc(self,pos):
    """Apply periodic boundary conditions

    Inputs
    ------
    pos : 1-d numpy array with length ndim
       a particle position which will be manipulated for PBCs

    """
    for idim in range(self.ndim):
      if pos[idim] < -self.half_boxlength[idim]: 
        pos[idim] += self.boxlength[idim]
      elif pos[idim] > self.half_boxlength[idim]: 
        pos[idim] -= self.boxlength[idim]
