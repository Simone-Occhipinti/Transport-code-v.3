# global variables
# space in cm
# energy in eV
import numpy as np

# GODIVA TEST

PARTICLE_TYPE = 'neutron'
#PARTICLE_TYPE = 'adjunction'
GEOMETRY_TYPE = 'sphere'
#GEOMETRY_TYPE = 'slab'
EMIN = 1E-5
EMAX = 2E7
EE = (EMIN,EMAX)
EREF = 2E6
Nstories = 1E5
SOURCE_POSITION = [0]
L0 = 0
LEnd = 8.7407
#LL = (L0,LEnd)
LL = (0.,)
WW = 6
Wmin = 1/WW
Wmax = WW
#Groups = [6.25E-1,5.53E3,8.21E5,1E6,EMAX]
Groups = np.logspace(np.log10(EMIN),np.log10(EMAX),int(1E3))

# only for eigenvalue calculationd
Nskip = 1E3
Nactive = 1E6
Kin = 1.