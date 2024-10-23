# global variables
# space in cm
# energy in eV

# GODIVA TEST

PARTICLE_TYPE = 'neutron'
GEOMETRY_TYPE = 'sphere'
EMIN = 1E-5
EMAX = 2E7
EE = (EMIN,EMAX)
EREF = 2E6
Nstories = 1E4
SOURCE_POSITION = [0]
L0 = 0
LEnd = 8.7407
LL = (L0,LEnd)
WW = 6
Wmin = 1/WW
Wmax = WW
Groups = [0.025,EMAX]

# only for eigenvalue calculationd
Nskip = 1E3
Nactive = 1E6
Kin = 1.