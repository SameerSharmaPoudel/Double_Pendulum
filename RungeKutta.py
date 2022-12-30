# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 12:42:28 2020

@author: Group 5
"""
from RKExplicit import ExplicitRK
from RKImplicit import ImplicitRK
def rungeKutta(RHS,t0,T,y0,dt,alpha,beta,gamma):
    num_ele_triu = countUpperTriangleNonZero(beta);
    if(num_ele_triu==0):
        return ExplicitRK(RHS,t0,T,y0,dt,alpha,beta,gamma);
    else:
        return ImplicitRK(RHS,t0,T,y0,dt,alpha,beta,gamma);
        
def countUpperTriangleNonZero(m):
    nonzero = 0;
    c = m.size;
    r = m[0].size;
    if(c==1):
        if m[0]==0:
            return 0
        else:
            return 1
    for i in range(0,c):
        for j in range(0+i,r):
            if(m[i][j]!=0):
                nonzero+=1;
    return nonzero;