# global variables
# space in cm
# energy in eV
import numpy as np

# GODIVA TEST

PARTICLE_TYPE = 'neutron'
GEOMETRY_TYPE = 'sphere'
EMIN = 1E-5
EMAX = 2E7
EE = (EMIN,EMAX)
EREF = 2E6
Nstories = 1E6
SOURCE_POSITION = [0]
L0 = 0
LEnd = 8.7407
LL = (0.,)
WW = 6
Wmin = 1/WW
Wmax = WW
#Groups = [0.025,EMAX]
Groups = np.logspace(np.log10(EMIN),np.log10(EMAX),int(1E3))

# only for eigenvalue calculationd
Nskip = 1E3
Nactive = 1E6
Kin = 1.