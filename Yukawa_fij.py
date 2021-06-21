#!/usr/bin/env python
# coding: utf-8

import numpy as np
import NEUTRINO2018

def ran_num():
    
    x = np.exp(np.random.uniform(np.log(1.*10**(-4)),np.log(10**(0))))
    
    return x

#H random matrix
def h():
    h = np.array([[ran_num(),ran_num(),ran_num()],
                  [ran_num(),ran_num(),ran_num()],
                  [ran_num(),ran_num(),ran_num()]])
    
    return h

#inverse of h matrix
def h1(h):
    
    h1 = np.linalg.inv(h)
    return h1

#matrix squared of neutrino masses
def D(m1,m2,m3):
    x = np.array([[np.sqrt(m1),0.,0.],
                  [0.,np.sqrt(m2),0.],[0.,0.,np.sqrt(m3)]])
    return x

def R_NH():
    x = np.array([[0.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,1.]])
    return x

def R_IH():
    x = np.array([[1.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,0.]])
    return x

def exp_val():
    
    k = NEUTRINO2018.NH()
    UPMNS = np.array([[k['U11'],k['U12'],k['U13']],
              [k['U21'],k['U22'],k['U23']],
              [k['U31'],k['U32'],k['U33']]])
    
    D_M = D(k['mnu1'],k['mnu2'],k['mnu3'])
    
    return UPMNS,D_M

#Yukawa without Lambda
def yuk_f():
    
    hij = h() 
    N = exp_val()
    fij = np.dot(np.dot(np.dot(np.dot(N[1],R_NH()),N[1]),np.transpose(N[0])),np.transpose(h1(hij)))
    return hij,fij,N,h1(hij)



