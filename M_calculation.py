# M calculation

import numpy as np
import os
from models import globalvariables as GV
from models import material_model as mat_m

base_dir = 'cross_sections_Janis'

U238_total = np.loadtxt(os.path.join(base_dir, 'U238', 'U238_total.csv'), delimiter=';', skiprows=3)
U238_scattering = np.loadtxt(os.path.join(base_dir, 'U238', 'U238_scattering.csv'), delimiter=';', skiprows=3)
U238_fission = np.loadtxt(os.path.join(base_dir, 'U238', 'U238_fission.csv'), delimiter=';', skiprows=3)
U238_nu = np.loadtxt(os.path.join(base_dir, 'U238', 'U238_nu.csv'), delimiter=';', skiprows=3)

U235_total = np.loadtxt(os.path.join(base_dir, 'U235', 'U235_total.csv'), delimiter=';', skiprows=3)
U235_scattering = np.loadtxt(os.path.join(base_dir, 'U235', 'U235_scattering.csv'), delimiter=';', skiprows=3)
U235_fission = np.loadtxt(os.path.join(base_dir, 'U235', 'U235_fission.csv'), delimiter=';', skiprows=3)
U235_nu = np.loadtxt(os.path.join(base_dir, 'U235', 'U235_nu.csv'), delimiter=';', skiprows=3)

U234_total = np.loadtxt(os.path.join(base_dir, 'U234', 'U234_total.csv'), delimiter=';', skiprows=3)
U234_scattering = np.loadtxt(os.path.join(base_dir, 'U234', 'U234_scattering.csv'), delimiter=';', skiprows=3)
U234_fission = np.loadtxt(os.path.join(base_dir, 'U234', 'U234_fission.csv'), delimiter=';', skiprows=3)
U234_nu = np.loadtxt(os.path.join(base_dir, 'U234', 'U234_nu.csv'), delimiter=';', skiprows=3)

uranium238 = mat_m.isotope(92,238,4.4984E21,U238_total[:,0],U238_total[:,1],U238_scattering[:,1],U238_fission[:,1],U238_nu[:,1])
uranium235 = mat_m.isotope(92,235,4.4994E22,U235_total[:,0],U235_total[:,1],U235_scattering[:,1],U235_fission[:,1],U235_nu[:,1])
uranium234 = mat_m.isotope(92,234,4.9184E20,U234_total[:,0],U234_total[:,1],U234_scattering[:,1],U234_fission[:,1],U234_nu[:,1])

MM = 1
S234 = np.trapz(uranium234.nu*uranium234.micro_xs_fission*uranium234.atomic_density,uranium234.energy)
f234 = lambda eout: uranium234.nu_avg(eout)*uranium234.macro_xs_fission(eout)/S234

S235 = np.trapz(uranium235.nu*uranium235.micro_xs_fission*uranium235.atomic_density,uranium235.energy)
f235 = lambda eout: uranium235.nu_avg(eout)*uranium235.macro_xs_fission(eout)/S235

S238 = np.trapz(uranium238.nu*uranium238.micro_xs_fission*uranium238.atomic_density,uranium238.energy)
f238 = lambda eout: uranium238.nu_avg(eout)*uranium238.macro_xs_fission(eout)/S238

gg = lambda eout: 1/(eout*(np.log(GV.EMAX)-np.log(GV.EMIN)))

out = 0

while out == 0:
    check = 0
    for ii in uranium235.energy:
        if MM*gg(ii) < f234(ii):
            check += 1
        if MM*gg(ii) < f235(ii):
            check += 1
        if MM*gg(ii) < f238(ii):
            check += 1
    if check == 0:
        out = MM
    else:
        MM += 1

print(MM)