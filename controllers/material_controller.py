# material controller

import numpy as np
from math import pi
import numpy.random as rnd
import models.material_model as mat
import models.physic_model as phy
import models.geometry_models as geo
import models.globalvariables as GV
import controllers.geometry_controller as geo_C

def find_energy_index(pp=float,mat=np.array):
    if pp >= mat[-1]:
        out = len(mat)-1
    elif pp <= mat[0]:
        out = 0
    else:
        out = np.where(pp<=mat)[0][0]
    return out

def find_position(pp=geo.point,dom=geo.domain):
    vett = []
    for ii in dom.materialposition:
        vett.append(ii[1])
    if pp.distance < GV.LEnd:
        index = np.where(pp.distance < np.array(vett))[0][0]
    else:
        index = len(dom.materialposition)-1
    return index

def sample_isotope_index(pp=phy.particle,mat=geo.domain):
    mat_index = find_position(pp.position,mat)
    nn = len(mat.materials[mat_index].composition)
    XS = np.zeros(nn)
    for ii in range(nn):
        XS[ii] += mat.materials[mat_index].composition[ii].macro_xs_total(pp.energy)/mat.materials[mat_index].macro_xs_total(pp.energy)
    cumXS = np.cumsum(XS)
    rho = rnd.rand()
    out = np.where(cumXS>rho)[0][0]
    return out

def sample_fission_isotope(pp=phy.particle,mat=mat.material):
    if mat.macro_xs_fission(pp.energy) > 0:
        nn = len(mat.composition)
        XS = np.zeros(nn)
        for ii in range(nn):
            XS[ii] += mat.composition[ii].macro_xs_fission(pp.energy)/mat.macro_xs_fission(pp.energy)
        cumXS = np.cumsum(XS)
        rho = rnd.rand()
        out = np.where(cumXS>rho)[0][0]
    else:
        out = None
    return out


