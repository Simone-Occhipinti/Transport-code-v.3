# statistic controller
import numpy as np
from math import pi
import models.globalvariables as GV
import models.statistic_model as stat
import models.physic_model as phy
import models.geometry_models as geo
import controllers.material_controller as mat_c
import numpy.random as rnd

def count_interaction(tt=stat.tally, pp=phy.particle, mat=geo.domain):
    if len(tt.spacerange) > 1:
        spaceindex = np.where(pp.position.distance <= tt.spacerange)[0][0]-1
    else:
        spaceindex = 0
    energyindex = np.where(pp.energy <= tt.energyrange)[0][0]-1
    mat_index = mat_c.find_position(pp.position,mat)
    if GV.simulation_type == 'placzek':
        xx = pp.weight
    else:
        xx = pp.weight/mat.materials[mat_index].macro_xs_total(pp.energy)
    tt.counter[energyindex][spaceindex] += xx

def wellford(tt=stat.tally):
    for ii in range(len(tt.mean)):
        delta = tt.counter[ii] - tt.mean[ii]
        tt.mean[ii] += delta/tt.iter
        delta2 = tt.counter[ii] - tt.mean[ii]
        tt.variance[ii] += delta*delta2
    tt.reset()

def normalization(tt=stat.tally):
    # normalizzo sull'estensione del gruppo energetico e sul volume del detector
    diffE = np.diff(tt.energyrange)
    for ii in range(len(diffE)):
        tt.mean[ii] *= 1/diffE[ii]
    if len(tt.spacerange)>1:
        for ii in range(len(tt.mean)):
            for jj in range(len(tt.mean[ii])):
                if GV.GEOMETRY_TYPE == 'sphere':
                    tt.mean[ii][jj] *= 1/(4/3*pi*(tt.spacerange[jj+1]**3-tt.spacerange[jj]**3))
                else:
                    tt.mean[ii][jj] *= 1/(tt.spacerange[jj+1]-tt.spacerange[jj])

def russian_roulette(pp=phy.particle):
    if pp.weight <= GV.Wmin:
        #rho = rnd.rand()
        rho = GV.rnd_counter.number()
        if rho <= pp.weight:
            pp.weight = 1
        else:
            pp.eof = 0
        
def splitting(pp=phy.particle,PS = list):
    if pp.weight >= GV.Wmax:
        N = pp.weight/GV.Wmax
        if N == int(N):
            N = int(N)
            for _ in range(N):
                PS.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/N))
            pp.weight *= 1/N
        else:
            D = N - int(N)
            #rho = rnd.rand()
            rho = GV.rnd_counter.number()
            if rho <= 1-D:
                for _ in range(int(N)-1):
                    PS.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/int(N)))
                pp.weight *= 1/int(N)
            else:
                for _ in range(int(N)):
                    PS.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/(int(N)+1)))
                pp.weight *= 1/(int(N)+1)

def restart_cycle(PS=list):
    PS = []

        
