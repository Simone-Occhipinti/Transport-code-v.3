# statistic models
import numpy as np
import models.globalvariables as GV
import models.physic_model as phy

class tally:
    def __init__(self,space_lim=tuple,space_n=int):
        # generazione della mesh spaziale
        if len(space_lim) > 1:
            self.spacerange = np.linspace(space_lim[0],space_lim[1],space_n)
            ll = len(self.spacerange)-1
        else:
            self.spacerange = np.array([0])
            ll = 1
        # generazione della mesh energetica
        self.energyrange = np.array(GV.Groups)
        self.mean = np.array([np.zeros(ll) for __ in range(len(self.energyrange))])
        self.variance = np.array([np.zeros(ll) for __ in range(len(self.energyrange))])
        self.iter = 0
    @property
    def avg(self):
        return self.mean
    @property
    def sigma(self):
        return np.sqrt(self.variance/(self.iter-1))
    
particle_squeue = []

        
        
        

