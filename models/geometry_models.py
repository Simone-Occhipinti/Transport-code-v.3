# Geometria e materiali

import numpy as np
from math import pi
import numpy.random as rnd
import models.globalvariables as GV
import models.material_model as mat

class point:
    def __init__(self, ics=float, ips=float, zeta=float):
        self.x = ics
        self.y = ips
        self.z = zeta
    @property
    def distance(self):
        if GV.GEOMETRY_TYPE == 'sphere':
            dd = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        else:
            dd = self.z
        return dd

class direction:
    def __init__(self, PHI=float, TE=float, degr=bool):
        if degr is True:
            PHI *= pi/180
            TE *= pi/180
        if TE >= 0 and TE <= pi:
            self.teta = TE
            self.phi = PHI
        elif PHI >= 0 and PHI <= pi:
            self.teta = PHI
            self.phi = TE
        else:
            self.teta = 'errore'
            self.phi = 'errore'
    def get_rnd_direction():
        phi = 2*pi*rnd.rand()
        teta = np.arccos(2*rnd.rand()-1)
        return direction(phi,teta,False)
    
class domain:
    def __init__(self, mate=list[tuple], space_ex=tuple, nspace=int, en_ex=tuple, nenergy=int, log=bool):
        # generazione disctretizzazione energia
        if log==True:
            self.energyrange = np.logspace(np.log10(en_ex[0]),np.log10(en_ex[1]),nenergy)
        else:
            self.energyrange = np.linspace(en_ex[0],en_ex[1],nenergy)
        # generazione della mesh spaziale
        if len(space_ex)>1:
            self.spacerange = np.linspace(space_ex[0],space_ex[1],nspace)
        else:
            self.spacerange = 1
        # generazione della distribuzione dei materiali
        self.materials = []
        self.materialposition = []
        for ii in range(len(mate)):
            self.materials.append(mate[ii][0])
            self.materialposition.append((mate[ii][1],mate[ii][2]))