# global variables

PARTICLE_TYPE = 'neutron'
GEOMETRY_TYPE = 'sphere'
EMIN = 1E-5
EMAX = 2E7
EE = (EMIN,EMAX)
SIMULATION_TYPE = 'source_given'
EREF = 2E6
Nstories = 1E6
SOURCE_POSITION = (0,0,0)
L0 = 0
LEnd = 300
LL = (L0,LEnd)
particle_squeue = []
WW = 8
Wmin = 1/WW
Wmax = WW
Groups = (EMIN, 1, 1E4, EMAX)