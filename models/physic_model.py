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
    def __init__(self,space_lim=tuple,space_n=int,initial_dist=list,intensity=float,type=str,pdf=None,range=None,pdf_log=None):
        self.intensity = intensity
        self.type = type
        self.shannonentropy = []
        if pdf != None:
            self.energydistribution = pdf
            if pdf_log == True:
                self.energyrange = np.logspace(np.log10(range[0]),np.log10(range[1]),len(pdf))
            else:
                self.energyrange = np.linspace(range[0],range[1],len(pdf))
        if len(space_lim)>1:
            self.spacerange = np.linspace(space_lim[0],space_lim[1],space_n)
            self.spacedistribution = np.zeros(len(self.spacerange),dtype=int)
        else:
            self.spacerange = np.array([float('inf')])
            self.spacedistribution = np.array([0])
        self.n_generated = np.zeros(len(self.spacerange),dtype=int)
        for ii in initial_dist:
            self.spacedistribution[ii] += 1
            self.n_generated[ii] += int(GV.Nstories/len(initial_dist))
    def get_position(self):
        indices = np.where(self.spacedistribution == 1)[0]
        rho = np.random.choice(indices)
        if len(self.spacerange)>1:
            out = self.spacerange[rho]
        else:
            out = 0
        return out
    def get_energy(self):
        if self.type == 'watt':
            out = watt(0.988,2.249)
        elif self.type == 'fixed':
            out = GV.EREF
        else:
            out = rejection(self.energydistribution,self.energyrange,1)
        if out > GV.EMAX:
            out = GV.EMAX
        elif out < GV.EMIN:
            out = GV.EMIN
        return out
    @property
    def tot_generated(self):
        return sum(self.n_generated)
    
    def s_entropy(self):
        HH = 0
        for ii in self.n_generated:
            if ii > 0:
                HH += -ii*np.log2(ii)
        self.shannonentropy.append(HH)
        
    def reset_source(self):
        self.spacedistribution = np.zeros(len(self.spacerange),dtype=int)
        self.n_generated = np.zeros(len(self.spacerange),dtype=int)

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