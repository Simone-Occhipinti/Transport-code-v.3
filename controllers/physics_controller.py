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

def generate_new_particle(ss=phy.source,ww=float,mat=mat.isotope,pp=None):
    if len(GV.LL) <= 1 or pp == None:
        rr = geo.point((0,0,0))
    else:
        rr = pp
    new_dir = geo.direction.get_rnd_direction()
    new_energy = ss.get_energy(mat)
    out = phy.particle(rr,new_dir,new_energy,ww)
    return out

def initialize_population(ss=phy.source,PS=list,mat=geo.domain):
    if all(len(sublist) == 0 for sublist in ss.spacedistribution):
        for _ in range(int(GV.Nstories)):
            if len(mat.materials)>1:
                mat_index = rnd.randint(0,len(mat.materials)-1)
            else:
                mat_index = 0
            if len(mat.materials[mat_index].composition)>1:
                is_index = rnd.randint(0,len(mat.materials[mat_index].composition)-1)
            else:
                is_index = 0
            PS.append(generate_new_particle(ss,1.,mat.materials[mat_index].composition[is_index]))
    else:
        ww = GV.Nstories/ss.tot_generated
        for ii in range(len(ss.spacedistribution)):
            if len(ss.spacedistribution[ii])>0:
                for jj in ss.spacedistribution[ii]:
                    for _ in range(jj.nn):
                        if ww <= GV.Wmin:
                            rho = rnd.rand()
                            if rho <= ww:
                                out = generate_new_particle(ss,1.,mat.materials[jj.mat_i].composition[jj.iso_i],jj.position)
                                PS.append(out)
                        else:
                            out = generate_new_particle(ss,ww,mat.materials[jj.mat_i].composition[jj.iso_i],jj.position)
                            PS.append(out)

def sample_free_flight(nn=phy.particle, mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    dir = geo_c.find_direction(nn)
    if dir > 0:
        mm = len(mat.materials)-mat_index
    else:
        mm = mat_index+1+len(mat.materials)
    bm = np.zeros(mm)
    rr = np.zeros(mm)
    if dir >= 0:
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
                ss = geo_c.distance_to_surface(nn,mat.materialposition[ii][0],dir)
            else:
                ss = geo_c.distance_to_surface(nn,mat.materialposition[ii][1],dir)
            rr[mat_index-ii] += ss-(rr[mat_index-ii-1] if ii<mat_index else 0)
        for ii in range(len(mat.materials)):
            if ii == 0:
                ss = geo_c.distance_to_surface(nn,mat.materialposition[ii][1],dir)
            else:
                ss = geo_c.distance_to_surface(nn,mat.materialposition[ii][0],dir)
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
    new_rr = ll+rr[region-2] if region > 1 else ll
    if GV.PARTICLE_TYPE == 'neutron':
        add_x = new_rr * np.sin(nn.direction.teta) * np.cos(nn.direction.phi)
        add_y = new_rr * np.sin(nn.direction.teta) * np.sin(nn.direction.phi)
        add_z = new_rr * np.cos(nn.direction.teta)
    else:
        add_x = new_rr * np.sin(-nn.direction.teta) * np.cos(nn.direction.phi + np.pi)
        add_y = new_rr * np.sin(-nn.direction.teta) * np.sin(nn.direction.phi + np.pi)
        add_z = new_rr * np.cos(-nn.direction.teta)
    return geo_c.sumpos(nn.position,geo.point((add_x,add_y,add_z)))

def alternative_free_flight(nn=phy.particle,mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    rho = rnd.rand()
    rr = -np.log(rho)/mat.materials[mat_index].macro_xs_total(nn.energy)
    add_x = rr * np.sin(nn.direction.teta) * np.cos(nn.direction.phi)
    add_y = rr * np.sin(nn.direction.teta) * np.sin(nn.direction.phi)
    add_z = rr * np.cos(nn.direction.teta)
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
        if low_i == up_i:
            sigma = np.array([mat.materials[mat_index].composition[is_index].micro_xs_scattering[up_i-1],mat.materials[mat_index].composition[is_index].micro_xs_scattering[up_i]])*mat.materials[mat_index].composition[is_index].atomic_density
            ee = np.array([mat.materials[mat_index].energy[low_i-1],mat.materials[mat_index].energy[low_i]])
        else:
            ee = mat.materials[mat_index].composition[is_index].energy[low_i-1:up_i]
            sigma = mat.materials[mat_index].composition[is_index].micro_xs_scattering[low_i-1:up_i]*mat.materials[mat_index].composition[is_index].atomic_density
        SS = np.trapz(sigma/ee,ee)
        ff = lambda eout: (mat.materials[mat_index].macro_xs_scattering(eout) / SS / eout) if not np.isclose(eout, 0) else 0
        new_energy = stat.rejection(ff,ee)
    return new_energy

def implicit_capture(nn=phy.particle, mat=geo.domain):
    mat_index = mat_c.find_position(nn.position,mat)
    is_index = mat_c.sample_isotope_index(nn,mat)
    new_weight = nn.weight*mat.materials[mat_index].composition[is_index].macro_xs_scattering(nn.energy)/mat.materials[mat_index].composition[is_index].macro_xs_total(nn.energy)
    if GV.PARTICLE_TYPE == 'adjunction':
        alfa = mat.materials[mat_index].composition[is_index].alpha
        low_i = mat_c.find_energy_index(nn.energy,mat.materials[mat_index].composition[is_index].energy)
        up_i = mat_c.find_energy_index(nn.energy/alfa,mat.materials[mat_index].composition[is_index].energy)
        if low_i == up_i:
            sigma = np.array([mat.materials[mat_index].composition[is_index].micro_xs_scattering[up_i-1],mat.materials[mat_index].composition[is_index].micro_xs_scattering[up_i]])**mat.materials[mat_index].composition[is_index].atomic_density
            ee = np.array([mat.materials[mat_index].energy[low_i-1],mat.materials[mat_index].energy[low_i]])
        else:
            ee = mat.materials[mat_index].composition[is_index].energy[low_i-1:up_i]
            sigma = mat.materials[mat_index].composition[is_index].micro_xs_scattering[low_i-1:up_i]*mat.materials[mat_index].composition[is_index].atomic_density
        SS = np.trapz(sigma/ee,ee)
        new_weight *= 1/(1-alfa)/mat.materials[mat_index].composition[is_index].macro_xs_scattering(nn.energy)*SS
    return new_weight

def add_fissionsite(pp=geo.point,isot=tuple,nn=float,ss=phy.source):
        if len(GV.LL) > 1:
            space_index = mat_c.find_energy_index(pp.distance,ss.spacerange)-1
        else:
            space_index = 0
        rho = rnd.rand()
        N = int(nn)
        if rho < nn-N:
            N+=1
        out = phy.fission_site(N,pp,isot)
        ss.spacedistribution[space_index].append(out)

def implicit_fission(nn=phy.particle,mat=geo.domain,ss=phy.source,KK=float):
    mat_index = mat_c.find_position(nn.position,mat)
    is_index = mat_c.sample_fission_isotope(nn,mat.materials[mat_index])
    if is_index == None:
        kn = 0
    else:
        SigmaF = mat.materials[mat_index].composition[is_index].macro_xs_fission(nn.energy)
        SigmaT = mat.materials[mat_index].composition[is_index].macro_xs_total(nn.energy)
        if GV.PARTICLE_TYPE == 'neutron':
            nu = mat.materials[mat_index].composition[is_index].nu_avg(nn.energy)
        else:
            SS = np.trapz(mat.materials[mat_index].composition[is_index].nu*mat.materials[mat_index].composition[is_index].micro_xs_fission*mat.materials[mat_index].composition[is_index].atomic_density,mat.materials[mat_index].composition[is_index].energy)
            chi = phy.watt_distribution(nn.energy)
            nu = chi/SigmaF*SS
        kn = nn.weight*(nu*SigmaF)/SigmaT
        #RR = kn/KK
        RR = kn
        add_fissionsite(nn.position,(mat_index,is_index),RR,ss)
        #nn.weight *= (1+nu*SigmaF/SigmaT)
    return kn/GV.Nstories
