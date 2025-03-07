# model physics

import numpy as np
import random
from math import pi
import models.globalvariables as GV
import numpy.random as rnd
import models.geometry_models as geo
import models.material_model as mat
import models.statistic_model as stat

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

class fission_site:
    def __init__(self,num=int,pos=geo.point,index=tuple):
        self.position = pos
        self.mat_i = index[0]
        self.iso_i = index[1]
        self.nn = num


class source:
    def __init__(self,nGen=int,space_lim=tuple,space_n=int,type=str,pdf=None,xx=None,pdf_log=None):
        self.type = type
        if pdf != None:
            self.energydistribution = pdf
            if pdf_log == True:
                self.energyrange = np.logspace(np.log10(xx[0]),np.log10(xx[1]),len(pdf))
            else:
                self.energyrange = np.linspace(xx[0],xx[1],len(pdf))
        if len(space_lim)>1:
            self.spacerange = np.linspace(space_lim[0],space_lim[1],space_n)
            self.spaceref = (self.spacerange[:-1] + self.spacerange[1:]) / 2
            self.spacedistribution = [[] for _ in range(len(self.spaceref))]
        else:
            self.spacerange = np.array([float('inf')])
            self.spaceref = np.array([0])
            self.spacedistribution = []

    def get_energy(self,mat=mat.isotope):
        if self.type == 'watt':
            if GV.PARTICLE_TYPE == 'neutron':
                out = watt()
            else:
                SS = np.trapz(mat.nu*mat.micro_xs_fission*mat.atomic_density,mat.energy)
                ff = lambda eout: mat.nu_avg(eout)*mat.macro_xs_fission(eout)/SS
                out = stat.rejection(ff,mat.energy,log=True)
        elif self.type == 'fixed':
            out = GV.EREF
        if out > GV.EMAX:
            out = GV.EMAX
        elif out < GV.EMIN:
            out = GV.EMIN
        return out
    
    @property
    def tot_generated(self):
        tot = 0
        for ii in range(len(self.spacedistribution)):
            for jj in self.spacedistribution[ii]:
                tot += jj.nn
        return tot
    
    @property
    def s_entropy(self):
        HH = 0
        for ii in range(len(self.spacedistribution)):
            pp = 0
            for jj in self.spacedistribution[ii]:
                pp += jj.nn/GV.Nstories
            HH += -pp*np.log2(pp) if pp > 0 else 0
        return HH
        
    def reset_source(self):
        self.spacedistribution = [[] for _ in range(len(self.spaceref))]


def watt_distribution(eout, aa=0.988, bb=2.249):
    xx = eout/1E6
    out = 0.4527*(np.exp(-xx/0.965))*np.sinh(np.sqrt(2.29*xx))
    return out

def watt(aa=0.988, bb=2.249):
        kk = 1 + bb/(8*aa)
        LL = (kk + np.sqrt(kk**2 - 1))/aa
        MM = aa*LL-1
        out = 0
        while out<=0:
            #xx = -np.log(rnd.rand())
            xx = -np.log(GV.rnd_counter.number())
            #yy = -np.log(rnd.rand())
            yy = -np.log(GV.rnd_counter.number())
            if ((yy-MM*(xx+1))**2)<=bb*LL*xx:
                out += LL*xx
        return out*1E6

