# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 13:05:48 2020

@author: Group 5
@return: res.tvec, res.yvec
"""
import numpy as np
from Simulationresults import Simulationresults
def ExplicitRK(RHS,t0,T,y0,dt,alpha,beta,gamma):
    m=len(beta)
    n=int((T-t0)/dt)+1
    nf=len(y0)
    tvec=np.arange(t0,T+dt,dt)
    #print(y0)
    yvec = np.zeros((n,nf))
    yvec[0,:]=y0
    k = np.zeros((m,nf))
    
    for j in range(0,n-1):
        prod=np.zeros((1,nf))
        for i in range(0,m):
            addn = np.zeros((1,nf))
            for l in range(0,m):
                if(l<i):
                    addn[0,:] = addn[0,:] + beta[i,l]*k[l,:]
            k[i,:] = RHS(tvec[j]+alpha[i]*dt,yvec[j,:]+dt*addn[0,:])
            prod[0,:] = prod[0,:] + gamma[i]*k[i,:]
        yvec[j+1,:]=yvec[j,:]+dt*prod[0,:]
    res = Simulationresults(tvec, yvec)
    return res
    