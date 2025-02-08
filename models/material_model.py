# material model
import numpy as np
from math import pi

class isotope:
    def __init__(self, atm_num=int, mas_num=int, dens=float, energy=np.array, tot=np.array, scatt=np.array, fiss=None, nu=None):
        self.atomic_number = atm_num
        self.mass_number = mas_num
        self.alpha = ((mas_num-1)/(mas_num+1))**2
        self.micro_xs_total = tot/1E24
        self.micro_xs_scattering = scatt/1E24
        self.micro_xs_absorption = self.micro_xs_total-self.micro_xs_scattering
        if fiss is not None:
            self.micro_xs_fission = fiss/1E24
            self.nu = nu
            self.fiss_check = True
        else:
            self.micro_xs_fission = np.zeros(len(energy))
            self.nu = np.zeros(len(energy))
            self.fiss_check = False
        self.micro_xs_capture = self.micro_xs_absorption-self.micro_xs_fission
        self.atomic_density = dens
        self.energy = energy
    
    def macro_xs_total(self, en=float):
        if en >= self.energy[-1]:
            out = self.micro_xs_total[-1]
        elif en <= self.energy[0]:
            out = self.micro_xs_total[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.micro_xs_total[index-1],self.micro_xs_total[index]])
        return out*self.atomic_density
    
    def macro_xs_scattering(self, en=float):
        if en >= self.energy[-1]:
            out = self.micro_xs_scattering[-1]
        elif en <= self.energy[0]:
            out = self.micro_xs_scattering[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.micro_xs_scattering[index-1],self.micro_xs_scattering[index]])
        return out*self.atomic_density
    
    def macro_xs_absorption(self, en=float):
        if en >= self.energy[-1]:
            out = self.micro_xs_absorption[-1]
        elif en <= self.energy[0]:
            out = self.micro_xs_absorption[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.micro_xs_absorption[index-1],self.micro_xs_absorption[index]])
        return out*self.atomic_density
    
    def macro_xs_fission(self, en=float):
        if en >= self.energy[-1]:
            out = self.micro_xs_fission[-1]
        elif en <= self.energy[0]:
            out = self.micro_xs_fission[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.micro_xs_fission[index-1],self.micro_xs_fission[index]])
        return out*self.atomic_density
    
    def nu_avg(self, en=float):
        if en >= self.energy[-1]:
            out = self.nu[-1]
        elif en <= self.energy[0]:
            out = self.nu[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.nu[index-1],self.nu[index]])
        return out

class material:
    def __init__(self, parts=list[isotope]):
        self.energy = parts[0].energy
        self.composition = parts
        nn = len(parts)
        self.macro_total = np.zeros(len(parts[0].micro_xs_total))
        self.macro_scattering = np.zeros(len(parts[0].micro_xs_total))
        self.macro_absorption = np.zeros(len(parts[0].micro_xs_total))
        self.macro_fission = np.zeros(len(parts[0].micro_xs_total))
        self.macro_capture = np.zeros(len(parts[0].micro_xs_total))
        self.nu = np.zeros(len(parts[0].micro_xs_total))
        for ii in range(nn):
            self.macro_total += parts[ii].micro_xs_total*parts[ii].atomic_density
            self.macro_scattering += parts[ii].micro_xs_scattering*parts[ii].atomic_density
            self.macro_absorption += parts[ii].micro_xs_absorption*parts[ii].atomic_density
            self.macro_fission += parts[ii].micro_xs_fission*parts[ii].atomic_density
            self.macro_capture += parts[ii].micro_xs_capture*parts[ii].atomic_density
            self.nu += parts[ii].nu*parts[ii].micro_xs_fission*parts[ii].atomic_density
        if np.any(self.macro_fission>0):
            self.nu *= 1/self.macro_fission
    def macro_xs_total(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_total[-1]
        elif en <= self.energy[0]:
            out = self.macro_total[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_total[index-1],self.macro_total[index]])
        return out
    def macro_xs_scattering(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_scattering[-1]
        elif en <= self.energy[0]:
            out = self.macro_scattering[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_scattering[index-1],self.macro_scattering[index]])
        return out
    def macro_xs_absorption(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_absorption[-1]
        elif en <= self.energy[0]:
            out = self.macro_absorption[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_absorption[index-1],self.macro_absorption[index]])
        return out
    def macro_xs_fission(self, en=float):
        if en >= self.energy[-1]:
            out = self.macro_fission[-1]
        elif en <= self.energy[0]:
            out = self.macro_fission[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.macro_fission[index-1],self.macro_fission[index]])
        return out
    def nu_avg(self, en=float):
        if en >= self.energy[-1]:
            out = self.nu[-1]
        elif en <= self.energy[0]:
            out = self.nu[0]
        else:
            index = np.where(self.energy>en)[0][0]
            out = np.interp(en, [self.energy[index-1],self.energy[index]], [self.nu[index-1],self.nu[index]])
        return out