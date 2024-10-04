# material controller
import numpy as np
from math import pi
import numpy.random as rnd
import models.material_model as mat
import models.physic_model as phy
import models.geometry_models as geo
import controllers.geometry_controller as geo_C

def find_energy_index(pp=float,mat=np.array):
    if pp >= mat[-1]:
        out = -1
    elif pp <= mat[0]:
        out = 0
    else:
        out = np.where(pp<=mat)[0][0]
    return out

def find_position(pp=geo.point,dom=geo.domain):
    vett = []
    for ii in dom.materialposition:
        vett.append(ii[1])
    index = np.where(pp.distance <= np.array(vett))[0][0]
    return index

def sample_isotope_index(pp=phy.particle,mat=geo.domain):
    mat_index = find_position(pp.position,mat)
    nn = len(mat.materials[mat_index].composition)
    XS = np.zeros(nn)
    en_index = find_energy_index(pp.energy,mat.materials[mat_index].energy)
    for ii in range(nn):
        XS[ii] += mat.materials[mat_index].composition[ii].micro_xs_scattering[en_index]*mat.materials[mat_index].composition[ii].atomic_density/mat.materials[mat_index].macro_scattering[en_index]
    cumXS = np.cumsum(XS)
    rho = rnd.rand()
    out = np.where(cumXS>rho)[0][0]
    return out

def sample_fission_isotope(pp=phy.particle,mat=mat.material):
    en_index = find_energy_index(pp.energy,mat.energy)
    if mat.macro_xs_fission(pp.energy) > 0:
        nn = len(mat.composition)
        XS = np.zeros(nn)
        for ii in range(nn):
            XS += mat.composition[ii].micro_xs_fission[en_index]*mat.composition[ii].atomic_density/mat.macro_fission[en_index]
        cumXS = np.cumsum(XS)
        rho = rnd.rand()
        out = np.where(cumXS>rho)[0][0]
    else:
        out = None
    return out


