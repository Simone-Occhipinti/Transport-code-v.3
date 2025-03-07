# Geometria e materiali

import numpy as np
from math import pi
import numpy.random as rnd
import models.globalvariables as GV
import models.material_model as mat

class point:
    def __init__(self, pos=tuple):
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
    @property
    def distance(self):
        if GV.GEOMETRY_TYPE == 'sphere':
            dd = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        else:
            dd = abs(self.z)
        return dd

class direction:
    def __init__(self, PHI=float, TE=float, degr=bool):
        if degr is True:
            PHI *= pi/180
            TE *= pi/180
        self.teta = TE
        self.phi = PHI
        
    def get_rnd_direction():
        #phi = 2*pi*rnd.rand()
        phi = 2*pi*GV.rnd_counter.number()
        #teta = np.arccos(2*rnd.rand()-1)
        teta = np.arccos(2*GV.rnd_counter.number()-1)
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
        if len(space_ex)>1:
            for ii in range(len(mate)):
                self.materials.append(mate[ii][0])
                self.materialposition.append((mate[ii][1],mate[ii][2]))
        else:
            self.materials.append(mate[0][0])
            self.materialposition.append((0,1E10))