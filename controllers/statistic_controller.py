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

def wellford(tt=stat.tally, pp=phy.particle, mat=geo.domain, placzeck=None):
    if len(tt.spacerange) > 1:
        spaceindex = np.where(pp.position.distance <= tt.spacerange)[0][0]-1
    else:
        spaceindex = 0
    energyindex = np.where(pp.energy <= tt.energyrange)[0][0]
    mat_index = mat_c.find_position(pp.position,mat)
    if placzeck == None:
        xx = pp.weight/mat.materials[mat_index].macro_xs_scattering(pp.energy)
    else:
        xx = pp.weight
    delta = xx - tt.mean[energyindex][spaceindex]
    tt.mean[energyindex][spaceindex] += delta/tt.iter
    tt.variance[energyindex][spaceindex] += delta*(xx - tt.mean[energyindex][spaceindex])

def normalization(tt=stat.tally):
    # normalizzo sull'estensione del gruppo energetico e sul volume del detector
    diffE = np.diff(np.array([GV.EMIN]+list(GV.Groups)))/GV.EREF
    for ii in range(len(diffE)):
        tt.mean[ii] *= 1/diffE[ii]
        tt.variance[ii] *= 1/diffE[ii]**2
    if len(tt.spacerange)>1:
        for ii in range(len(tt.mean)):
            for jj in range(len(tt.mean[ii])):
                if GV.GEOMETRY_TYPE == 'sphere':
                    tt.mean[ii][jj] *= 1/(4/3*pi*(tt.spacerange[jj+1]**2-tt.spacerange[jj]**2))
                    tt.variance[ii][jj] *= 1/(4/3*pi*(tt.spacerange[jj+1]**2-tt.spacerange[jj]**2))**2
                else:
                    tt.mean[ii][jj] *= 1/(tt.spacerange[jj+1]-tt.spacerange[jj])
                    tt.variance[ii][jj] *= 1/(tt.spacerange[jj+1]-tt.spacerange[jj])**2

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
                stat.particle_squeue.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/N))
            pp.weight *= 1/N
        else:
            D = N - int(N)
            if rnd.rand() <= 1-D:
                for ii in range(int(N)-1):
                    stat.particle_squeue.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/int(N)))
                pp.weight *= 1/int(N)
            else:
                for ii in range(int(N)):
                    stat.particle_squeue.append(phy.particle(pp.position,pp.direction,pp.energy,pp.weight/(int(N)+1)))
                pp.weight *= 1/(int(N)+1)

def restart_cycle():
    stat.particle_squeue = []

        
