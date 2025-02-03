# global variables
# space in cm
# energy in eV
import numpy as np

# GODIVA TEST

#PARTICLE_TYPE = 'neutron'
PARTICLE_TYPE = 'adjunction'
GEOMETRY_TYPE = 'sphere'
#GEOMETRY_TYPE = 'slab'
EREF = 2E6
EMIN = 1E-5
#EMIN = EREF*np.exp(-4)
#EMIN = EREF
#EMAX = EREF
EMAX = 2E7
#EMAX = 4*EMIN
EE = (EMIN,EMAX)
Nstories = 1E1
SOURCE_POSITION = [(0.,0.,0.)]
L0 = 0
LEnd = 8.7407  #cm
#LL = (0.,)
Wmin = 0.25
Wmax = 20
#Groups = np.logspace(np.log10(EMIN),np.log10(EMAX),int(200))

# only for eigenvalue calculationd
simulation_type = 'eigenvalue'
#simulation_type = 'placzek'
Nskip = 1E2
Nactive = 1E6
Kin = 1.
Groups = [EMIN, 6.25E-1,5.53E3,8.21E5,1E6,EMAX]
LL = (L0,LEnd)