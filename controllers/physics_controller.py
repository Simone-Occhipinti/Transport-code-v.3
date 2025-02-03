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

def generate_new_particle(ss=phy.source,ww=float,mat=geo.domain,pp=None):
    if pp == None:
        rr = ss.get_position()
    else:
        rr = pp
    mat_index = mat_c.find_position(pp,mat)
    new_dir = geo.direction.get_rnd_direction()
    new_energy = ss.get_energy(mat.materials[mat_index])
    out = phy.particle(rr,new_dir,new_energy,ww)
    return out

def generate_population(ss=phy.source,ww=float,PS=list,mat=geo.domain):
    for ii in range(len(ss.spacedistribution)):
        for _ in range(ss.n_generated[ii]):
            PS.append(generate_new_particle(ss,ww,mat,ss.spacedistribution[ii]))

def sample_free_flight(nn=phy.particle, mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    dir = geo_c.find_direction(nn)
    if dir > 0:
        mm = len(mat.materials)-mat_index
    else:
        mm = mat_index+1+len(mat.materials)
    bm = np.zeros(mm)
    rr = np.zeros(mm)
    if dir > 0:
        for ii in range(mm):
            if ii == 0:
                ss = geo_c.distance_to_surface(nn,mat.materialposition[ii+mat_index][1],dir)
            else:
                ss = geo_c.distance_to_surface(nn,mat.materialposition[ii+mat_index][0],dir)
            rr[ii] += ss-(rr[ii-1] if ii>0 else 0)
        for ii in range(mm):
            bm[ii] = np.sum([mat.materials[mat_index + jj].macro_xs_total(nn.energy) * rr[jj] for jj in range(ii + 1)])
    else:
        for ii in range(mat_index,-1,-1):
            if ii == mat_index:
                ss = np.max([0,geo_c.distance_to_surface(nn,mat.materialposition[ii][0],dir)])
            else:
                ss = np.max([0,geo_c.distance_to_surface(nn,mat.materialposition[ii][1],dir)])
            rr[mat_index-ii] += ss-(rr[mat_index-ii-1] if ii<mat_index else 0)
        for ii in range(len(mat.materials)):
            if ii == 0:
                ss = np.max([0,geo_c.distance_to_surface(nn,mat.materialposition[ii][1],-dir)])
            else:
                ss = np.max([0,geo_c.distance_to_surface(nn,mat.materialposition[ii][0],-dir)])
            rr[mat_index+ii+1] += ss-rr[mat_index+ii]
        for ii in range(mat_index,-1,-1):
            bm[ii] = np.sum([mat.materials[ii-jj].macro_xs_total(nn.energy)*rr[jj] for jj in range(ii+1)])
        for ii in range(len(mat.materials)):
            bm[mat_index+ii+1] = np.sum([mat.materials[jj].macro_xs_total(nn.energy)*rr[mat_index+jj] for jj in range(ii+1)])
    rho = rnd.rand()
    eta = -np.log(rho)
    for region in range(1,mm+1):
        b_minus1 = bm[region - 2] if region > 1 else 0
        if b_minus1 < eta <= bm[region-1]:
            break
    if dir > 0:
        index = region-1
    else:
        index = mat_index-region+1 if region-1 <= mat_index else region-2-mat_index
    ll = (eta-b_minus1)/mat.materials[index].macro_xs_total(nn.energy)
    new_rr = ll+b_minus1
    if GV.PARTICLE_TYPE == 'neutron':
        add_x = new_rr * np.sin(nn.direction.teta) * np.cos(nn.direction.phi)
        add_y = new_rr * np.sin(nn.direction.teta) * np.sin(nn.direction.phi)
        add_z = new_rr * np.cos(nn.direction.teta)
    else:
        add_x = new_rr * np.sin(-nn.direction.teta) * np.cos(nn.direction.phi + np.pi)
        add_y = new_rr * np.sin(-nn.direction.teta) * np.sin(nn.direction.phi + np.pi)
        add_z = new_rr * np.cos(-nn.direction.teta)
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
        ee = mat.materials[mat_index].composition[is_index].energy[low_i:up_i+1]
        sigma = mat.materials[mat_index].macro_scattering[low_i:up_i+1]
        SS = np.trapz(sigma/ee,ee)
        ff = lambda eout: (mat.materials[mat_index].macro_xs_scattering(eout) / SS / eout) if not np.isclose(eout, 0) else 0
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
        ss.spacedistribution.append(pp)
        rho = rnd.rand()
        N = int(nn)
        if rho < nn-N:
            ss.n_generated.append(N+1)
        else:
            ss.n_generated.append(N)

def implicit_fission(nn=phy.particle,mat=geo.domain,ss=phy.source):
    mat_index = mat_c.find_position(nn.position,mat)
    is_index = mat_c.sample_fission_isotope(nn,mat.materials[mat_index])
    if is_index == None:
        kn = 0
    else:
        SigmaF = mat.materials[mat_index].macro_xs_fission(nn.energy)
        SigmaT = mat.materials[mat_index].macro_xs_total(nn.energy)
        if GV.PARTICLE_TYPE == 'neutron':
            nu = mat.materials[mat_index].nu_avg(nn.energy)
        else:
            SS = np.trapz(mat.materials[mat_index].nu*mat.materials[mat_index].macro_fission,mat.materials[mat_index].energy)
            chi = phy.watt_distribution(nn.energy)
            nu = chi/SigmaF*SS
        kn = nn.weight*(nu*SigmaF)/SigmaT
        add_fissionsite(nn.position,kn,ss)
    return kn/GV.Nstories
