# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:45:47 2020

@author: Group 5
"""
import numpy as np
def fixed_point(PHI,x0,TOL,Nmax):
    xn=x0
    for i in range(1,Nmax):
        xn1 = PHI(xn)
        if(np.linalg.norm(xn1-xn)<TOL):
            xn = xn1
            break
        xn = xn1
    return xn
