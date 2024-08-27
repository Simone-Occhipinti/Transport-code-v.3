# statistic controller
import numpy as np
from math import pi
import models.globalvariables as GV
import models.statistic_model as stat
import models.physic_model as phy
import models.geometry_models as geo
import controllers.physics_controller as phy_c
import controllers.material_controller as mat_c
import numpy.random as rnd

def new_iter(tt=stat.tally):
    tt.iter += 1

def new_interaction(tt=stat.tally,pp=phy.particle,mat=geo.domain):
    if len(tt.spacerange) > 1:
        spaceindex = np.where(pp.position.distance <= tt.spacerange)[0][0]-1
    else:
        spaceindex = 0
    energyindex = mat_c.find_energy_index(pp.energy,tt.energyrange) - 1
    mat_index = mat_c.find_position(pp.position,mat)
    tt.estimator[spaceindex][energyindex] += pp.weight/mat.materials[mat_index].macro_xs_scattering(pp.energy)
    tt.mom2[spaceindex][energyindex] += (pp.weight/mat.materials[mat_index].macro_xs_scattering(pp.energy))**2

def generate_sample_average(tt=stat.tally,type=str,n_groups=tuple):
    out1 = [np.zeros(len(tt.spacerange)-1) for __ in range(len(n_groups)-1)]
    out2 = [np.zeros(len(tt.spacerange)-1) for __ in range(len(n_groups)-1)]
    low_i = int(0)
    up_i = int(0)
    # raggruppo i valori per energie all'interno dello stesso gruppo
    for ii in range(len(n_groups)-1):
        up_i = mat_c.find_energy_index(n_groups[ii+1],tt.energyrange)
        for jj in range(len(tt.spacerange)-1):
            out1[ii][jj] += np.sum(tt.estimator[jj][low_i:up_i])/tt.iter
            out2[ii][jj] += np.sum(tt.mom2[jj][low_i:up_i])/tt.iter
        low_i = up_i
    out1 = np.array(out1)
    out2 = np.array(out2)
    diffE = np.diff(np.array(n_groups)/GV.EREF)
    # normalizzo sull'estensione del gruppo energetico e sul volume del detector
    for ii in range(len(diffE)):
        out1[ii] *= 1/diffE[ii]
        out2[ii] *= 1/diffE[ii]**2
    if len(tt.spacerange)>1:
        for ii in out1:
            for jj in range(len(ii)):
                if GV.GEOMETRY_TYPE == 'sphere':
                    ii[jj] *= 1/(4/3*pi*(tt.spacerange[jj+1]**2-tt.spacerange[jj]**2))
                else:
                    ii[jj] *= 1/(tt.spacerange[jj+1]-tt.spacerange[jj])
    var = out2 - out1**2
    std = np.sqrt(var/tt.iter)
    return (out1,std)

def russian_roulette(pp=phy.particle):
    if pp.weight <= GV.Wmin:
        rho = rnd.rand()
        if rho <= GV.Wmin:
            pp.weight *= GV.WW
        else:
            pp.eof = 0
        
def splitting(pp=phy.particle):
    if pp.weight >= GV.Wmax:
        N = pp.weight/GV.Wmax
        if N == int(N):
            for ii in range(N-1):
                GV.particle_squeue.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/N))
            pp.weight *= 1/N
        else:
            D = N - int(N)
            if rnd.rand() <= 1-D:
                for ii in range(int(N)-1):
                    GV.particle_squeue.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/int(N)))
                pp.weight *= 1/int(N)
            else:
                for ii in range(int(N)):
                    GV.particle_squeue.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/(int(N)+1)))
                pp.weight *= 1/(int(N)+1)
        

        
