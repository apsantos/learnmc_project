class inputs(object):
  def __init__(self):
    """Simulation inputs
    Can either be assigned in python or on command line
    """
    self.T = 2.
    self.mu = 0.1
    self.NVT = True
    self.muVT = False

    self.boxlength = [10., 10., 10.]
    self.ndim = len(self.boxlength)
    self.rc = min(self.boxlength)/2.

    self.nmol = int(108) # initial N for muVT
    if self.NVT:
      self.max_nmol = self.nmol # needed for muVT
    elif self.muVT:
      self.max_nmol = 10000 # needed for muVT

    if self.nmol > self.max_nmol:
        print('Max number of molecules is: %d' % self.max_nmol)
        sys.exit(main())

    self.nmoves = 20000
    self.nmovesequil = 10000
    self.nsample = 50

    # Move mixture
    if self.NVT:
      self.frac_displacement = 1.0
      self.frac_insdel = 0.0
    elif self.muVT:
      self.frac_displacement = 0.6
      self.frac_insdel = 0.4

    self.sigma = 1.
    self.epsilon = 1.

    self.make_movie = False
    self.fname = 'mc'

  def readcommandline(self,parser):
    """Read command line information to override the defaults

    Parameters
    ----------
    parser : an ArgParser.parser
        Carries command-line information

    """
    if parser.parse_args().V != 'read':
      V = float(parser.parse_args().V)
      v3 = V**(1./3.)
      self.boxlength  = [v3, v3, v3]
      self.rc = min(self.boxlength)/2.
    if parser.parse_args().N != 'read':
      self.nmol = int(parser.parse_args().N)
      self.max_nmol = self.nmol # needed for muVT
    if parser.parse_args().temperature != 'read':
      self.T = float(parser.parse_args().temperature)
    if parser.parse_args().mu != 'read':
      self.mu = float(parser.parse_args().mu)
      self.max_nmol = 10000
      self.NVT = False
      self.muVT = True
    else:
      self.muVT = True
      self.NVT = False
    if parser.parse_args().outfile != 'read':
      self.fname = parser.parse_args().outfile

