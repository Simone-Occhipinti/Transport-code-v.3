# statistic models
import numpy as np
import models.globalvariables as GV
import models.physic_model as phy

class tally:
    def __init__(self,space_lim=tuple,space_n=int,energy_lim=tuple,energy_n=int,log=bool):
        # generazione della mesh spaziale
        if len(space_lim) > 1:
            self.spacerange = np.linspace(space_lim[0],space_lim[1],space_n)
            self.deltaspace = np.diff(self.spacerange)
        else:
            self.spacerange = 0
            self.deltaspace = 0
        # generazione della mesh energetica
        self.energyrange = np.logspace(np.log10(energy_lim[0]),np.log10(energy_lim[1]),energy_n)
        self.energydelta = np.diff(self.energyrange)
        self.estimator = np.array([np.zeros(len(self.energyrange)-1) for __ in range(len(self.spacerange)-1)])
        self.mom2 = np.array([np.zeros(len(self.energyrange)-1) for __ in range(len(self.spacerange)-1)])
        self.iter = 0
        
        
        

