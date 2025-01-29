# Controller Physics
import numpy as np
import numpy.random as rnd

import controllers.geometry_controller as geo_c
import controllers.material_controller as mat_c
import models.physic_model as phy
import models.material_model as mat
import models.geometry_models as geo
import models.globalvariables as GV
import models.statistic_model as stat

def generate_new_particle(ss=phy.source,ww=float,pp=None):
    new_energy = ss.get_energy()
    if pp == None:
        rr = ss.get_position()
    else:
        rr = pp
    new_position = geo_c.get_point_from_surface(rr)
    if new_position.distance > GV.LEnd:
        new_position.x += -1E-4
        new_position.y += -1E-4
        new_position.z += -1E-4
    new_dir = geo.direction.get_rnd_direction()
    out = phy.particle(new_position,new_dir,new_energy,ww)
    return out

def generate_population(ss=phy.source,ww=float,PS=list):
    indices = np.where(ss.spacedistribution==1)[0]
    for ii in indices:
        for jj in range(ss.n_generated[ii]):
            PS.append(generate_new_particle(ss,ww,ss.spaceref[ii]))

def sample_free_flight(nn=phy.particle, mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    dir = geo_c.find_direction(nn)
    if dir > 0:
        mm = len(mat.materials)-mat_index
    else:
        mm = mat_index+1
    bm = np.zeros(mm)
    rr = np.zeros(mm)
    if dir > 0:
        for ii in range(mm):
            if ii == 0:
                ss = geo_c.distance_par2surf(nn,mat.materialposition[ii+mat_index][1])
            else:
                ss = geo_c.distance_par2surf(nn,mat.materialposition[ii+mat_index][0])
            rr[ii] += ss-rr[ii-1]
        for ii in range(mm):
            bm[ii] = np.sum([mat.materials[mat_index+jj].macro_xs_total(nn.energy)*rr[jj] for jj in range(ii+1)])
    else:
        for ii in range(mat_index,-1,-1):
            if ii == mat_index:
                ss = geo_c.distance_par2surf(nn,mat.materialposition[ii][0])
            else:
                ss = geo_c.distance_par2surf(nn,mat.materialposition[ii][1])
            rr[mat_index-ii] += ss-rr[mat_index-ii-1]
        for ii in range(mat_index,-1,-1):
            bm[ii] = np.sum([mat.materials[ii-jj].macro_xs_total(nn.energy)*rr[jj] for jj in range(ii+1)])
    rho = rnd.rand()
    eta = -np.log(rho)
    for region in range(1,mm+1):
        b_minus1 = bm[region-2] if ii>1 else 0
        if b_minus1 < eta <= bm[region-1]:
            break
    if dir > 0:
        index = region-1
    else:
        index = mat_index-region+1
    ll = (eta-b_minus1)/mat.materials[index].macro_xs_total(nn.energy)
    new_rr = ll+b_minus1
    if GV.PARTICLE_TYPE == 'neutron':
        add_x = new_rr*np.sin(nn.direction.teta)*np.cos(nn.direction.phi)
        add_y = new_rr*np.sin(nn.direction.teta)*np.sin(nn.direction.phi)
        add_z = new_rr*np.cos(nn.direction.teta)
    else:
        add_x = new_rr*np.sin(-nn.direction.teta)*np.cos(-nn.direction.phi)
        add_y = new_rr*np.sin(-nn.direction.teta)*np.sin(-nn.direction.phi)
        add_z = new_rr*np.cos(-nn.direction.teta)
    return geo_c.sumpos(nn.position,geo.point((add_x,add_y,add_z)))

def sample_energy_stepf(nn=phy.particle,mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    is_index = mat_c.sample_isotope_index(nn,mat)
    alfa = mat.materials[mat_index].composition[is_index].alpha
    if GV.PARTICLE_TYPE == 'neutron':
        new_energy = nn.energy*(alfa+(1-alfa)*rnd.rand())
    else:
        low_i = mat_c.find_energy_index(nn.energy,mat.materials[mat_index].composition[is_index].energy)
        up_i = mat_c.find_energy_index(nn.energy/alfa,mat.materials[mat_index].composition[is_index].energy)
        ee = mat.materials[mat_index].composition[is_index].energy[low_i:up_i]
        sigma = mat.materials[mat_index].macro_scattering[low_i:up_i]
        SS = np.trapz(sigma/ee,ee)
        ff = lambda eout: mat.materials[mat_index].macro_xs_scattering(eout)/SS/eout
        new_energy = stat.rejection(ff,ee)
    return new_energy

def new_weight(nn=phy.particle, mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    new_weight = nn.weight*mat.materials[mat_index].macro_xs_scattering(nn.energy)/mat.materials[mat_index].macro_xs_total(nn.energy)
    if GV.PARTICLE_TYPE == 'adjunction':
        is_index = mat_c.sample_isotope_index(nn,mat)
        alfa = mat.materials[mat_index].composition[is_index].alpha
        low_i = mat_c.find_energy_index(nn.energy,mat.materials[mat_index].composition[is_index].energy)
        up_i = mat_c.find_energy_index(nn.energy/alfa,mat.materials[mat_index].composition[is_index].energy)
        SS = np.trapz(mat.materials[mat_index].macro_scattering[low_i:up_i]/mat.materials[mat_index].composition[is_index].energy[low_i:up_i],mat.materials[mat_index].composition[is_index].energy[low_i:up_i])
        new_weight *= 1/(1-alfa)/mat.materials[mat_index].macro_xs_scattering(nn.energy)*SS
    return new_weight

def add_fissionsite(pp=geo.point,nn=float,ss=phy.source):
        space_index = np.where(ss.spacerange>=pp.distance)[0][0]
        ss.spacedistribution[space_index] = 1
        rho = rnd.rand()
        N = int(nn)
        if rho < nn-N:
            ss.n_generated[space_index] += N+1
        else:
            ss.n_generated[space_index] += N

def implicit_fission(nn=phy.particle,mat=geo.domain,kk=float,ss=phy.source):
    mat_index = mat_c.find_position(nn.position,mat)
    is_index = mat_c.sample_fission_isotope(nn,mat.materials[mat_index])
    if is_index == None:
        kn = 0
    else:
        en_index = mat_c.find_energy_index(nn.energy,mat.materials[mat_index].composition[is_index].energy)
        #SigmaF = mat.materials[mat_index].composition[is_index].micro_xs_fission[en_index]*mat.materials[mat_index].composition[is_index].atomic_density
        SigmaF = mat.materials[mat_index].macro_xs_fission(nn.energy)
        SigmaT = mat.materials[mat_index].macro_xs_total(nn.energy)
        RR = nn.weight*(mat.materials[mat_index].composition[is_index].nu*SigmaF)/SigmaT/kk
        kn = RR*kk/GV.Nstories
        add_fissionsite(nn.position,RR,ss)
    return kn
