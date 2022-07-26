# python code for monte carlo simulations
import sys, argparse
import numpy as np
from box import box
from configs import configfile
from configs import initialize
from configs import movie
from energy import energy
from input import inputs
from moves import moves
from tally import tally

def run(parser):
#def run(parser=argparse.ArgumentParser()):
  inp = inputs()
  inp.readcommandline(parser)
  
  bx = box()
  bx.setedges(inp.boxlength)
  
  en = energy(inp.nmol, bx, inp.T)
  en.setparameters(inp.sigma,inp.epsilon,inp.rc)
  
  if inp.muVT:
    init = initialize(inp.max_nmol)
  else:
    init = initialize(inp.nmol)
  pos, exist = init.random(inp.nmol, bx)
  
  # Set up move
  mv = moves(inp.nmol,bx,en,inp.T,inp.mu)
  mv.setmovemix(inp.frac_displacement,inp.frac_insdel)
  
  if inp.make_movie:
    mov = movie(inp.sigma,bx)
    mov.update(pos,exist)
  
  config = configfile(inp.ndim, inp.fname)
  io = tally(inp.fname)
  
  for imove in range(inp.nmoves):
    mv.move(pos, exist)
    if (imove % inp.nsample) == 0 and (imove > inp.nmovesequil):
      en.calcpressure(pos, mv.nmol, mv.vol)
      en.calcfullenergy(pos, mv.nmol)
      config.writeframe(imove, pos, mv.nmol)
      io.update(imove, en, mv)
  
  if inp.make_movie:
    mov.update(pos, exist)
  io.average()
  # io.pltseries()
  if inp.muVT:
    io.outputave(inp.T, inp.fname,inp.mu)
  else:
    io.outputave(inp.T, inp.fname)
  mv.outputefficiency(inp.fname)

def main(argv=None):
  parser = argparse.ArgumentParser(description='Simple MC code')
  parser.add_argument("-V", type=str, default='read',
                 help='volume, assumes')
  parser.add_argument("-N", type=str, default='read',
                 help='Number of particles')
  parser.add_argument('-T', "--temperature", type=str, default='read',
                 help='Temperature')
  parser.add_argument("--mu", type=str, default='read',
                 help='chemical potential')
  parser.add_argument("--outfile", type=str, default='read',
                 help='output file name root')

  run(parser)

if __name__ == '__main__':
    sys.exit(main())
