# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:51:06 2020

@author: Group 5
"""

from RungeKutta import rungeKutta
import numpy as np
def ExplicitEuler(RHS,t0,T,y0,dt):
    alpha = np.array([0.0])
    beta = np.array([[0.0]])
    gamma = np.array([1.0])
    return rungeKutta(RHS, t0, T, y0, dt, alpha, beta, gamma)

def ImplicitEuler(RHS,t0,T,y0,dt):
    alpha = np.array([1.0])
    beta = np.array([[1.0]])
    gamma = np.array([1.0])
    return rungeKutta(RHS, t0, T, y0, dt, alpha, beta, gamma)

def RK4(RHS,t0,T,y0,dt):
    alpha = np.array([0, 0.5, 0.5, 1])
    beta = np.array([[0, 0, 0, 0],[0.5, 0, 0, 0],[0, 0.5, 0, 0],[0, 0, 1, 0]])
    gamma = np.array([1/6, 1/3, 1/3, 1/6])
    return rungeKutta(RHS, t0, T, y0, dt, alpha, beta, gamma)

def Midpoint(RHS,t0,T,y0,dt):
    alpha = np.array([0.5])
    beta = np.array([[0.5]])
    gamma = np.array([1.0])
    return rungeKutta(RHS, t0, T, y0, dt, alpha, beta, gamma)