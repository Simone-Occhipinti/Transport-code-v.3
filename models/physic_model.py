# model physics

import numpy as np
from math import pi
import models.globalvariables as GV
import numpy.random as rnd
import models.geometry_models as geo

class particle:
    def __init__(self, pos=geo.point, dir=geo.direction, ee=float, ww=float):
        self.position = pos
        self.direction = dir
        self.energy = ee
        self.neutron_mass = 1.675E-27
        @property
        def velocity(self):
            out = np.sqrt(2*self.energy/self.neutron_mass)
            return out 
        self.weight = ww
        self.eof = 1

class source:
    def __init__(self, points=list[geo.point],intensity=float, type=str, pdf=None, range=None, pdf_log=None):
        self.intensity = intensity
        self.type = type
        if pdf != None:
            self.distribution = pdf
            if pdf_log == True:
                self.range = np.logspace(np.log10(range[0]),np.log10(range[1]),len(pdf))
            else:
                self.range = np.linspace(range[0],range[1],len(pdf))
        self.position = points
    def get_position(self):
        if len(self.position)==1:
            return self.position[0]
        else:
            rho = np.ceil(rnd.rand()*len(self.position))
            return self.position[rho]
    def get_energy(self):
        if self.type == 'watt':
            return watt(0.988,2.249)
        else:
            return rejection(self.distribution,self.range,1)

def watt(aa, bb):
        kk = 1 + bb/(8*aa)
        LL = (kk + np.sqrt(kk**2 - 1))/aa
        MM = aa*LL-1
        out = 0
        while out<=0:
            xx = -np.log(rnd.rand())
            yy = -np.log(rnd.rand())
            if ((yy-MM*(xx+1))**2)<=bb*LL*xx:
                out += LL*xx
        return out*1E6

def rejection(fun=np.array, var=np.array, nn=int):
    if len(fun) != len(var):
        return 'length_error'
    if np.trapz(fun,var) != 1:
        return 'error_not_normalized'
    aa, bb = var[0], var[-1]
    max = np.max(fun)
    out = []
    for ii in range(nn):
        rho = 0
        while rho == 0:
            pp1, pp2 = aa + (bb-aa)*rnd.rand(), max*rnd.rand()
            index = np.where(pp1>=var)
            if index == 1:
                index+=1
            refernce = np.interp(pp1, [var[index-1], var[index]], [fun[index-1], fun[index]])
            if pp2 <= refernce:
                out.append(pp2)
    return out