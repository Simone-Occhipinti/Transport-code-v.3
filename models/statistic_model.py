# statistic models
import numpy as np
import numpy.random as rnd
import models.globalvariables as GV
import models.physic_model as phy

class tally:
    def __init__(self,space_lim=tuple,space_n=int):
        # generazione della mesh spaziale
        if len(space_lim) > 1:
            self.spacerange = np.linspace(space_lim[0],space_lim[1],space_n)
            self.spaceref = (self.spacerange[:-1] + self.spacerange[1:]) / 2
            ll = len(self.spaceref)
        else:
            self.spacerange = np.array([0])
            self.spaceref = np.array([0])
            ll = 1
        # generazione della mesh energetica
        self.energyrange = np.array(GV.Groups)
        self.energyref = (self.energyrange[:-1] + self.energyrange[1:]) / 2
        self.mean = np.array([np.zeros(ll) for __ in range(len(self.energyref))])
        self.variance = np.array([np.zeros(ll) for __ in range(len(self.energyref))])
        self.iter = 0
        self.counter = np.array([np.zeros(ll) for __ in range(len(self.energyref))])
    
    def reset(self):
        if len(self.spacerange)>1:
            ll = len(self.spaceref)
        else:
            ll = 1
        self.counter = np.array([np.zeros(ll) for __ in range(len(self.energyref))])

    @property
    def avg(self):
        return self.mean
    @property
    def sigma(self):
        return np.sqrt((self.variance/self.iter))
    @property
    def sigma_smavg(self):
        aa = self.sigma
        return aa/np.sqrt(self.iter)
    @property
    def RSD(self):
        aa = self.sigma_smavg
        return aa/self.mean

def rejection(ff, vett=np.array):
    yy = []
    for ii in vett:
        yy.append(ff(ii))
    yy = np.array(yy)
    MM = np.max(yy)
    out = 0
    while out == 0:
        ics = vett[0] + (vett[-1]-vett[0])*rnd.rand()
        ips = MM*rnd.rand()
        if ips <= ff(ics):
            out = ics
    return out

        
        
        

