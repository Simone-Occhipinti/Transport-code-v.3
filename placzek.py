import numpy as np
import matplotlib.pyplot as plt
from models import globalvariables as GV


E0 = GV.EREF

def placzek1(aa):
    alfa = ((aa-1)/(aa+1))**2
    e0 = 2E6
    xx = np.linspace(alfa*e0,e0,100)
    yy = ((e0/xx)**(alfa/(1-alfa)))/(xx*(1-alfa))
    return [xx,yy]

def placzek2(aa):
    alfa = ((aa-1)/(aa+1))**2
    e0 = 2E6
    xx = np.linspace(alfa*alfa*e0,e0*alfa,100)
    yy = (((e0/xx)**(alfa/(1-alfa)))/(xx*(1-alfa)**2))*((1-alfa)*(1-(alfa**(1/(1-alfa)))) - (alfa**(alfa/(1-alfa)))*np.log(alfa*e0/xx))
    return [xx,yy]

def placzek3(aa):
    alfa = ((aa-1)/(aa+1))**2
    e0 = 2e6
    csi = 1+((alfa*np.log(alfa))/(1-alfa))
    xx = np.linspace(1E-5,e0*alfa*alfa,100)
    yy = 1/csi/xx
    return [xx,yy]

def placzek(AA):
    part1 = placzek1(AA)
    part2 = placzek2(AA)
    part3 = placzek3(AA)

    xlab = np.hstack((part3[0],part2[0],part1[0]))
    ylab = np.hstack((part3[1],part2[1],part1[1]))

    return [xlab,ylab]


def adj_placzek1(aa,fsigma):
    alfa = ((aa-1)/(aa+1))**2
    xx = np.linspace(E0,E0/alfa,200)
    sigma = []
    for ii in xx:
        sigma.append(fsigma(ii))
    sigma = np.array(sigma)
    sigma0 = fsigma(E0)
    #sigma = 1
    #yy = (1/(1-alfa))*(xx/E0)**(1/(1-alfa)) 
    yy = sigma*(1/(E0*sigma0*(1-alfa)))*(xx/E0)**(alfa/(1-alfa))
    return [xx,yy]

def adj_placzek2(aa,fsigma):
    alfa = ((aa-1)/(aa+1))**2
    xx = np.linspace(E0/alfa,E0/(alfa**2),200)
    sigma = []
    for ii in xx:
        sigma.append(fsigma(ii))
    sigma = np.array(sigma)
    sigma0 = fsigma(E0)
    #yy = (1/(1-alfa)**2)*((xx/E0)**(1/(1-alfa)))*(((1-alfa)*(1-alfa**(1/(1-alfa))))-((alfa**(1/(1-alfa)))*(np.log(alfa*xx/E0))))
    #sigma = 1
    yy = sigma*(1/(E0*sigma0*(1-alfa)**2))*((xx/E0)**(alfa/(1-alfa)))*(((1-alfa)*(1-alfa**(1/(1-alfa))))-((alfa**(1/(1-alfa)))*np.log(alfa*xx/E0)))
    return [xx,yy]

def adj_placzek3(aa,fsigma):
    alfa = ((aa-1)/(aa+1))**2
    xx = np.linspace(E0/alfa**2,E0/(alfa**3),200)
    sigma = []
    for ii in xx:
        sigma.append(fsigma(ii))
    sigma = np.array(sigma)
    sigma0 = fsigma(E0)
    #sigma = 1
    yy = sigma*(1/(E0*sigma0*(1-alfa)**3))*((xx/E0)**(alfa/(1-alfa)))*((((1-alfa)**2)*(1-alfa**(1/(1-alfa))))-((1-alfa)*(alfa**(1/(1-alfa)))*np.log(alfa*xx/E0))+((alfa**(2/(1-alfa)))*(((1-alfa)*np.log(xx/E0*alfa**2))+(((np.log(xx/E0*alfa**2))**2)/2))))
    return [xx,yy]

def adj_placzek(AA,sigma):
    part1 = adj_placzek1(AA,sigma)
    part2 = adj_placzek2(AA,sigma)
    part3 = adj_placzek3(AA,sigma)

    xlab = np.hstack((part1[0],part2[0],part3[0]))
    ylab = np.hstack((part1[1],part2[1],part3[1]))

    return [xlab,ylab]


# letargy


def Lplaczek1(aa):
    alfa = ((aa-1)/(aa+1))**2
    eps = np.log(1/alfa)
    xx = np.linspace(0,eps,100)
    yy = 1/(1-alfa)*np.exp((alfa/(1-alfa))*xx)
    return [xx,yy]

def Lplaczek2(aa):
    alfa = ((aa-1)/(aa+1))**2
    eps = np.log(1/alfa)
    xx = np.linspace(eps,2*eps,100)
    yy = (np.exp((alfa/(1-alfa))*xx))/(1-alfa) * (1 - ((1-((eps-xx)/(1-alfa)))*np.exp(-eps/(1-alfa))))
    return [xx,yy]

def Lplaczek3(aa):
    alfa = ((aa-1)/(aa+1))**2
    eps = np.log(1/alfa)
    csi = 1 + alfa/(1-alfa)*np.log(alfa)
    xx = np.linspace(2*eps,25,100)
    yy = np.array([1/csi for _ in range(len(xx))])
    return [xx,yy]

def Lplaczek(AA):
    part1 = Lplaczek1(AA)
    part2 = Lplaczek2(AA)
    part3 = Lplaczek3(AA)

    xlab = np.hstack((part1[0],part2[0],part3[0]))
    ylab = np.hstack((part1[1],part2[1],part3[1]))

    return [xlab,ylab]
