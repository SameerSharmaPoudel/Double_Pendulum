# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 13:05:48 2020

@author: Group 5
@return: res.tvec, res.yvec
"""
import numpy as np
from Simulationresults import Simulationresults
from NonlinearSolvers import fixed_point

def ImplicitRK(RHS,t0,T,y0,dt,alpha,beta,gamma):
    n = y0.size
    m = gamma.size
    steps = int((T-t0)/dt)+1
    tvec = np.zeros((steps,1))
    yvec = np.zeros((steps,n))
    tvec[0] = t0
    yvec[0,:] = y0
    for i in range(1,steps):
        ti = t0+i*dt
        kgsum = np.zeros((1, n))
        k0 = np.zeros(m)
        for mi in range(0,m):
            gammai = gamma[mi]
            alphai = alpha[mi]
            betai = beta[:,mi];
            betaisum = np.ones((1,m))*betai;
            def PHI(k):
                return RHS(ti+alphai*dt, y0 + (dt * k * betaisum)[0,:])
            ki = fixed_point(PHI, k0, 1E-12, 30)
            kgsum = kgsum + gammai * ki;
        y1=(y0 + dt * kgsum)[0,:];
        tvec[i] = ti;
        yvec[i,:] = y1;
        y0 = y1;
    res = Simulationresults(tvec, yvec)
    return res

    