#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess
import time
import sys

#PMNS matrix 2018 to 3\sigma https://arxiv.org/abs/1708.01186
def nuexpvalues():   

    dcOut={}
    #phases of the PMNS matrix and the R 
    phases1 = np.random.uniform(0.,0.0*np.pi,3) # WARNING! They are in zero
    delta = phases1[0]
    eta1 = phases1[1]
    eta2 = phases1[2]

    #light neutrino masses (up 3 sigma range) NH (NO)
    #mnu1 = 10**((np.log10(2.5e-3)-np.log10(1e-9))*np.random.uniform(0,1)+np.log10(1e-9))*1e-9 
    #Zero in this framework
    mnu1 = 1.0e-20
    mnu2 = np.sqrt(np.random.uniform(7.05e-5,8.14e-5)*1.0e-18+mnu1**2)
    mnu3 = np.sqrt(np.random.uniform(2.41e-3,2.60e-3)*1.0e-18+mnu1**2)

    #mixing angles (up 3 sigma range) NH
    t12 = np.arcsin(np.sqrt(np.random.uniform(0.273,0.379)))
    t23 = np.arcsin(np.sqrt(np.random.uniform(0.445,0.599)))
    t13 = np.arcsin(np.sqrt(np.random.uniform(0.0196,0.0241)))

    #Building PMNS matrix NH
    UM12 = np.array([ [np.cos(t12),np.sin(t12),0.], [-np.sin(t12),np.cos(t12),0.], [0.,0.,1.0] ])
    UM13 = np.array([ [np.cos(t13),0.,np.sin(t13)], [0.,1.0,0.], [-np.sin(t13),0.,np.cos(t13)] ])
    UM23 = np.array([ [1.0,0.,0.], [0.,np.cos(t23),np.sin(t23)], [0.,-np.sin(t23),np.cos(t23)] ])
    Uphases = np.array([ [np.exp(eta1*1j),0.,0.], [0.,np.exp(eta2*1j),0.], [0.,0.,1.0] ])
    U=np.dot(UM23,np.dot(UM13,np.dot(UM12,Uphases)))

    dcOut['mnu1']= mnu1
    dcOut['mnu2']= mnu2
    dcOut['mnu3']= mnu3
    dcOut['U11']= np.real(U[0,0])*1.
    dcOut['U12']= np.real(U[0,1])*1.
    dcOut['U13']= np.real(U[0,2])*1.
    dcOut['U21']= np.real(U[1,0])*1.
    dcOut['U22']= np.real(U[1,1])*1.
    dcOut['U23']= np.real(U[1,2])*1.
    dcOut['U31']= np.real(U[2,0])*1.
    dcOut['U32']= np.real(U[2,1])*1.
    dcOut['U33']= np.real(U[2,2])*1.

    return dcOut
