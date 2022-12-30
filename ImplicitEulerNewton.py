# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:52:37 2020

@author: Group 5
"""
import numpy as np
import sympy as sp
from Simulationresults import Simulationresults

def ImplicitEulerNewton(dt,t_end,cfg,y0):
    
    '''Input - length, masses, initial momenta and positions of both masses, acceleration due to gravity, 
               time step size and total time
        Returns Position and momenta at every time step'''
    l1,l2,m1,m2,g = cfg.l1, cfg.l2, cfg.m1, cfg.m2, cfg.g
    p1_initial,p2_initial,theta1_initial,theta2_initial=y0[2],y0[3],y0[0],y0[1]

    theta1_new,theta2_new,p1_new,p2_new = [],[],[],[]
    theta1_new.append(theta1_initial)
    theta2_new.append(theta2_initial)
    p1_new.append(p1_initial)
    p2_new.append(p2_initial)
    
    N = int(t_end/dt)+1
    tvec = np.zeros(N)
    
    t1n,t2n,p1n,p2n,t1o,t2o,p1o,p2o = sp.symbols('t1n,t2n,p1n,p2n,t1o,t2o,p1o,p2o')
    X = [t1n,t2n,p1n,p2n,t1o,t2o,p1o,p2o]
    
    c = sp.cos(t1n - t2n)
    s = sp.sin(t1n - t2n)
    
    r1 = (t1n - t1o - dt*(((p1n*l2) - (p2n*l1*c))/ ((l1**2)*l2*(m1+m2*(s**2)))))
    r2 = (t2n - t2o - dt*(((p2n*(m1+m2)*l1)-(m2*p1n*l2*c))/(m2*l1*(l2**2)*(m1+m2*(s**2)))))  

    temp1 = 1/((l1**2)*(l2**2)*((((c**2)*m2)-m1-m2)**2))
    temp2 = g*sp.sin(t1n)*(l1**3)*(l2**2)*(m2**2)*(m1+m2)*(c**4)
    temp3 = 2*l1*m2*l2*(c**2)*((-0.5*p1n*p2n*s)+(g*sp.sin(t1n)*(l1**2)*l2*((m1+m2)**2)))
    temp4 = s*c*(((p2n**2)*(l1**2)*(m1+m2))+((l2**2)*(p1n**2)*m2))
    temp5 = l1*l2*(m1+m2)*((p1n*p2n*s)+(g*sp.sin(t1n)*(l1**2)*l2*((m1+m2)**2)))
    
    temp6 = g*(l1**2)*(l2**3)*(m2**3)*sp.sin(t2n)*(c**4)
    temp7 = 2*l1*m2*l2*(c**2)* ((0.5*p1n*p2n*s)+(g*sp.sin(t2n)*(l2**2)*l1*(m1+m2)*m2))
    temp8 = l1*l2*(m1+m2)*((-p1n*p2n*s)+(g*sp.sin(t2n)*(l2**2)*l1*(m1+m2)*m2))
    
    r3 = (p1n - p1o - (temp1*(-temp2+temp3+temp4-temp5))*dt)
    r4 = (p2n - p2o - (temp1*(-temp6+temp7-temp4-temp8))*dt)

    Rr = sp.Matrix([r1,r2,r3,r4])
    
    r11 = sp.diff(r1,t1n)
    r12 = sp.diff(r1,t2n)
    r13 = sp.diff(r1,p1n)
    r14 = sp.diff(r1,p2n)
    r21 = sp.diff(r2,t1n)
    r22 = sp.diff(r2,t2n)
    r23 = sp.diff(r2,p1n)
    r24 = sp.diff(r2,p2n)
    r31 = sp.diff(r3,t1n)
    r32 = sp.diff(r3,t2n)
    r33 = sp.diff(r3,p1n)
    r34 = sp.diff(r3,p2n)
    r41 = sp.diff(r4,t1n)
    r42 = sp.diff(r4,t2n)
    r43 = sp.diff(r4,p1n)
    r44 = sp.diff(r4,p2n)

    Gr = sp.Matrix([[r11,r12,r13,r14],[r21,r22,r23,r24],[r31,r32,r33,r34],[r41,r42,r43,r44]])
    
    R = sp.lambdify(X, Rr,'numpy')
    G = sp.lambdify(X, Gr,'numpy')
    
    
    t1o,t2o,p1o,p2o,t1n,t2n,p1n,p2n = theta1_initial,theta2_initial,p1_initial,p2_initial,0,0,0,0
    
    new = [t1n,t2n,p1n,p2n]
    for i in range(N):
        tvec[i] = dt*i
        normres=1
        tol = 10**(-7)
        counter = 0
        maxiter = 5000
        #print("\n")
        while((normres>tol)and(counter<=maxiter)):
            temp = G(t1n,t2n,p1n,p2n,t1o,t2o,p1o,p2o)
            temp = np.linalg.inv(temp)
            dot = np.dot(temp,R(t1n,t2n,p1n,p2n,t1o,t2o,p1o,p2o))
            for j in range(4):
                new[j]=float(new[j]-dot[j])
            t1n,t2n,p1n,p2n = new[0],new[1],new[2],new[3]
            normres=np.linalg.norm(R(t1n,t2n,p1n,p2n,t1o,t2o,p1o,p2o))
            #print("The norm of the residuum vector is {}".format(normres))
            counter+=1
        t1o,t2o,p1o,p2o = new[0],new[1],new[2],new[3]
        if(i<N-1):
            theta1_new.append(new[0])
            theta2_new.append(new[1])
            p1_new.append(new[2])
            p2_new.append(new[3])
        if counter>maxiter:
            print("Solution did not converge for {}-th step".format(i+1))
    yvec = np.array([theta1_new, theta2_new, p1_new, p2_new]).transpose()
    res = Simulationresults(tvec, yvec)
    return res