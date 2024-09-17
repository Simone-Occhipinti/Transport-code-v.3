# geometry controller
import numpy as np
from math import pi
import models.globalvariables as GV
import models.geometry_models as geo
import models.physic_model as phy

def sumpos(pp=geo.point,qq=geo.point):
    xx = pp.x + qq.x
    yy = pp.y + qq.y
    zz = pp.z + qq.z
    out = geo.point((xx,yy,zz))
    return out

def difpos(pp=geo.point,qq=geo.point):
    xx = pp.x - qq.x
    yy = pp.y - qq.y
    zz = pp.z - qq.z
    out = geo.point((xx,yy,zz))
    return out

def find_direction(pp=phy.particle):
    delta_x = np.sin(pp.direction.teta)*np.cos(pp.direction.phi)
    delta_y = np.sin(pp.direction.teta)*np.sin(pp.direction.phi)
    delta_z = np.cos(pp.direction.teta)
    delta = geo.point((delta_x,delta_y,delta_z))
    prova = sumpos(pp.position,delta)
    if pp.position.distance <= prova.distance:
        return 1
    else:
        return -1

def distance_par2surf(pp=phy.particle,ll=float):
    r_punto = pp.position.z*np.cos(pp.direction.teta) + pp.position.y*np.sin(pp.direction.teta)*np.sin(pp.direction.phi) + pp.position.x*np.sin(pp.direction.teta)*np.cos(pp.direction.phi)
    deltaD = ll**2 - pp.position.distance**2
    Delta = r_punto**2 + deltaD**2
    return -r_punto+Delta

def is_outofbound(pp=phy.particle,type=str):
    if type == 'space':
        if pp.position.distance > GV.LEnd or pp.position.distance < GV.L0:
            return True
        else:
            return False
    else:
        if pp.energy > GV.EMAX or pp.energy < GV.EMIN:
            return True
        else:
            return False