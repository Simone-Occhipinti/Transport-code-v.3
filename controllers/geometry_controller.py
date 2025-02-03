# geometry controller
import numpy as np
import numpy.random as rnd
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

def distance_to_surface(pp=phy.particle, ll=float, direction=int):

    xx, yy, zz= pp.position.x, pp.position.y, pp.position.z
    theta, phi = pp.direction.teta, pp.direction.phi

    # Convert spherical angles to unit vector
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    
    if GV.GEOMETRY_TYPE == 'sphere':
        AA = dx**2 + dy**2 + dz**2
        BB = 2*(dx*xx + dy*yy + dz*zz)
        CC = xx**2 + yy**2 + zz**2 - ll**2

        Delta = BB**2 - 4*AA*CC
        if Delta >= 0:
            t1 = (-BB - np.sqrt(Delta))/(2*AA)
            t2 = (-BB + np.sqrt(Delta))/(2*AA)

            x1, y1, z1 = xx+dx*t1, yy+dy*t1, zz+dz*t1
            x2, y2, z2 = xx+dx*t2, yy+dy*t2, zz+dz*t2

            D1 = np.sqrt(((xx-x1)**2 + (yy-y1)**2 + (zz-z1)**2))
            D2 = np.sqrt(((xx-x2)**2 + (yy-y2)**2 + (zz-z2)**2))

            out = min([D1,D2]) if direction > 0 else max([D1,D2])
        else:
            out = 0

    elif GV.GEOMETRY_TYPE == 'slab':
        out = abs((ll-zz)/dz)
        
    return out

def is_outofbound(pp=phy.particle,type=str):
    if type == 'space':
        if len(GV.LL)>1:
            if pp.position.distance > GV.LEnd or pp.position.distance < GV.L0:
                return True
            else:
                return False
        else:
            return False
    else:
        if pp.energy > GV.EMAX or pp.energy < GV.EMIN:
            return True
        else:
            return False

def get_point_from_surface(rr=float):
    if GV.GEOMETRY_TYPE == 'sphere':
        dir = geo.direction.get_rnd_direction()
        xx = rr*np.sin(dir.teta)*np.cos(dir.phi)
        yy = rr*np.sin(dir.teta)*np.sin(dir.phi)
        zz = rr*np.cos(dir.teta)
    else:
        xx = 0
        yy = 0
        zz = rr
    return geo.point((xx,yy,zz))