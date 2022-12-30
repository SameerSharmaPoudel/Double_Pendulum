# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:04:33 2020

@author: Group 5
"""

import numpy as np
import sympy as sp
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from IPython.display import HTML
import logging

logger = logging.getLogger()

def theta1dot(t, x, cfg):
    t1 = x[0]
    t2 = x[1]
    p1 = x[2]
    p2 = x[3]
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    
    c = np.cos(t1 - t2)
    s = np.sin(t1 - t2)
    theta1dot = ((p1*l2) - (p2*l1*c))/((l1**2)*l2*(m1+m2*(s**2)))
    return theta1dot

def theta2dot(t, x, cfg):
    t1 = x[0]
    t2 = x[1]
    p1 = x[2]
    p2 = x[3]
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    
    c = np.cos(t1 - t2)
    s = np.sin(t1 - t2)
    theta2dot = ((p2*(m1+m2)*l1)-(m2*p1*l2*c))/ (m2*l1*(l2**2)*(m1+m2*(s**2)))  
    return theta2dot

def p1dot(t, x, cfg):
    t1 = x[0]
    t2 = x[1]
    p1 = x[2]
    p2 = x[3]
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    g = cfg.g
    
    c = np.cos(t1 - t2)
    s = np.sin(t1 - t2)
    
    temp1 = 1/((l1**2)*(l2**2)*((((c**2)*m2)-m1-m2)**2))
    
    temp2 = g*np.sin(t1)*(l1**3)*(l2**2)*(m2**2)*(m1+m2)*(c**4)
    temp3 = 2*l1*m2*l2*(c**2)*((-0.5*p1*p2*s)+g*np.sin(t1)*(l1**2)*l2*((m1+m2)**2))
    temp4 = s*c*(((p2**2)*(l1**2)*(m1+m2))+((l2**2)*(p1**2)*m2))
    temp5 = l1*l2*(m1+m2)*((p1*p2*s)+(g*np.sin(t1)*(l1**2)*l2*((m1+m2)**2)))
        
    p1dot = temp1*(-temp2+temp3+temp4-temp5)
    return p1dot

def p2dot(t, x, cfg):
    t1 = x[0]
    t2 = x[1]
    p1 = x[2]
    p2 = x[3]
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    g = cfg.g
    
    c = np.cos(t1 - t2)
    s = np.sin(t1 - t2)
    
    temp1 = 1/((l1**2)*(l2**2)*((((c**2)*m2)-m1-m2)**2))
    temp4 = s*c*(((p2**2)*(l1**2)*(m1+m2))+((l2**2)*(p1**2)*m2))
    temp6 = g*(l1**2)*(l2**3)*(m2**3)*np.sin(t2)*(c**4)
    temp7 = 2*l1*m2*l2*(c**2)*((0.5*p1*p2*s)+(g*np.sin(t2)*(l2**2)*l1*(m1+m2)*m2))
    temp8 = l1*l2*(m1+m2)*((-p1*p2*s)+(g*np.sin(t2)*(l2**2)*l1*(m1+m2)*m2))
    p2dot = temp1*(-temp6+temp7-temp4-temp8)
    return p2dot

def DoublePendulum (t,x, cfg):
    p1dt = p1dot(t,x, cfg)
    p2dt = p2dot(t,x, cfg)
    theta1dt = theta1dot(t,x, cfg)
    theta2dt = theta2dot(t,x, cfg)
    return np.array([theta1dt, theta2dt, p1dt, p2dt])

def DoublePendulumPrime (t,x,cfg):    
    t1n,t2n,p1n,p2n = sp.symbols('t1n,t2n,p1n,p2n')
    X = [t1n,t2n,p1n,p2n]
    
    c = sp.cos(t1n - t2n)
    s = sp.sin(t1n - t2n)
    
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    g = cfg.g
    
    tdot1 = ((p1n*l2) - (p2n*l1*c))/ ((l1**2)*l2*(m1+m2*(s**2)))
    tdot2 = ((p2n*(m1+m2)*l1)-(m2*p1n*l2*c))/(m2*l1*(l2**2)*(m1+m2*(s**2)))

    temp1 = 1/((l1**2)*(l2**2)*((((c**2)*m2)-m1-m2)**2))
    temp2 = g*sp.sin(t1n)*(l1**3)*(l2**2)*(m2**2)*(m1+m2)*(c**4)
    temp3 = 2*l1*m2*l2*(c**2)*((-0.5*p1n*p2n*s)+(g*sp.sin(t1n)*(l1**2)*l2*((m1+m2)**2)))
    temp4 = s*c*(((p2n**2)*(l1**2)*(m1+m2))+((l2**2)*(p1n**2)*m2))
    temp5 = l1*l2*(m1+m2)*((p1n*p2n*s)+(g*sp.sin(t1n)*(l1**2)*l2*((m1+m2)**2)))
    
    temp6 = g*(l1**2)*(l2**3)*(m2**3)*sp.sin(t2n)*(c**4)
    temp7 = 2*l1*m2*l2*(c**2)* ((0.5*p1n*p2n*s)+(g*sp.sin(t2n)*(l2**2)*l1*(m1+m2)*m2))
    temp8 = l1*l2*(m1+m2)*((-p1n*p2n*s)+(g*sp.sin(t2n)*(l2**2)*l1*(m1+m2)*m2))
    
    pdot1 = (temp1*(-temp2+temp3+temp4-temp5))
    pdot2 = (temp1*(-temp6+temp7-temp4-temp8))
    
    d11 = sp.diff(tdot1,t1n)
    d12 = sp.diff(tdot1,t2n)
    d13 = sp.diff(tdot1,p1n)
    d14 = sp.diff(tdot1,p2n)
    d21 = sp.diff(tdot2,t1n)
    d22 = sp.diff(tdot2,t2n)
    d23 = sp.diff(tdot2,p1n)
    d24 = sp.diff(tdot2,p2n)
    d31 = sp.diff(pdot1,t1n)
    d32 = sp.diff(pdot1,t2n)
    d33 = sp.diff(pdot1,p1n)
    d34 = sp.diff(pdot1,p2n)
    d41 = sp.diff(pdot2,t1n)
    d42 = sp.diff(pdot2,t2n)
    d43 = sp.diff(pdot2,p1n)
    d44 = sp.diff(pdot2,p2n)

    Gr = sp.Matrix([[d11,d12,d13,d14],[d21,d22,d23,d24],[d31,d32,d33,d34],[d41,d42,d43,d44]])
    
    G = sp.lambdify(X, Gr,'numpy')
    t1n,t2n,p1n,p2n = x[0],x[1],x[2],x[3]
    
    return G(t1n,t2n,p1n,p2n)

def get_positions(res, cfg):
    X1 = cfg.l1*np.sin(res.yvec[:,0])
    Y1 = -cfg.l1*np.cos(res.yvec[:,0])
    X2 = X1+cfg.l2*np.sin(res.yvec[:,1])
    Y2 = Y1-cfg.l2*np.cos(res.yvec[:,1])
    return X1, Y1, X2, Y2

def plot_positions(res, cfg):
    fig, ax = plt.subplots()
    maxsize=cfg.l1+cfg.l2
    ax.set(xlim=(-maxsize, maxsize), ylim=(-maxsize, maxsize))
    ax.set_aspect('equal')
    
    
    X1,Y1,X2,Y2 = get_positions(res, cfg)
    
    ax.scatter(X1, Y1,s=1,label='Mass 1')
    ax.scatter(X2, Y2,s=1,label='Mass 2')
    
    ax.legend()
    ax.set_xlabel('X, m')
    ax.set_ylabel('Y, m')
    ax.set_title('Position curve through time')
    fig.show()
    
    fig, ax = plt.subplots()
    
    ax.plot(res.tvec, res.yvec[:,0],label='Mass 1')
    ax.plot(res.tvec, res.yvec[:,1],label='Mass 2')
    ax.legend()
    ax.set_ylabel('Position, rad')
    ax.set_xlabel('Time, s')
    ax.set_title('Position')
    fig.show()

def plot_velocities(res, cfg):
    fig, ax = plt.subplots()
    
    
    t1dot, t2dot = get_velocities(res, cfg)
    
    ax.plot(res.tvec, t1dot,label='Mass 1')
    ax.plot(res.tvec, t2dot,label='Mass 2')
    ax.legend()
    
    ax.set_xlabel('Time, s')
    ax.set_ylabel('Angular velocity, rad/s')
    ax.set_title('Velocities')
    fig.show()
    
def plot_accelerations(res, cfg):
    fig, ax = plt.subplots()
    
    
    t1acc, t2acc = get_accelerations(res, cfg)
    
    ax.plot(res.tvec, t1acc,label='Mass 1')
    ax.plot(res.tvec, t2acc,label='Mass 2')
    
    ax.legend()
    ax.set_xlabel('Time, s')
    ax.set_ylabel('Angular acceleration, rad/s^2')
    ax.set_title('Accelerations')
    fig.show()

def plot_energy(res, cfg):
    fig, ax = plt.subplots()
    #ax.legend()
        
    E = get_energy(res, cfg)
           
    #maxsize=np.amax(E)
    #minsize=np.amin(E)
    #minmultip=0.9
    #maxmultip=1.1
    #if(maxsize<0):
    #    maxmultip = 0.9
    #if(minsize<0):
    #    minmultip = 1.1
        
    #ax.set(ylim=(minsize*minmultip, maxsize*maxmultip))
    
    #ax.plot(res.tvec, E)
    #ax.set_xlabel('Time, s')
    #ax.set_ylabel('Total energy')
    #ax.set_title('Total energy')
    #fig.show()  
    plt.plot(res.tvec, E)
    plt.ylabel('Total energy')
    plt.xlabel('Time, s')
    plt.title('Total energy')
    plt.show()

def get_velocities(res, cfg):
    
    '''Input - Position and Momenta at all time instances, length and mass
       Returns velocities at every time instance'''
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    
    t1 = res.yvec[:,0]
    t2 = res.yvec[:,1]
    p1 = res.yvec[:,2]
    p2 = res.yvec[:,3]
    t1dot,t2dot = [0]*len(t1),[0]*len(t2)
    for i in range(len(t1)):
        c = np.cos(t1[i] - t2[i])
        s = np.sin(t1[i] - t2[i])
        t1dot[i] = ((p1[i]*l2) - (p2[i]*l1*c))/ ((l1**2)*l2*(m1+m2*(s**2)))
        t2dot[i] = ((p2[i]*(m1+m2)*l1)-(m2*p1[i]*l2*c))/ (m2*l1*(l2**2)*(m1+m2*(s**2))) 
        
    return t1dot, t2dot

def get_accelerations(res, cfg):
    
    '''Input - Position and Velocities at all time instances, length and mass and g
       Returns accelerations at every time instance'''
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    g = cfg.g
    
    t1 = res.yvec[:,0]
    t2 = res.yvec[:,1]
    t1dot,t2dot = get_velocities(res, cfg)
    
    t1_acc,t2_acc = [0]*len(t1),[0]*len(t2)
    for i in range(len(t1)):
        c = np.cos(t1[i] - t2[i])
        s = np.sin(t1[i] - t2[i])
        t1_acc[i] = (m2*g*np.sin(t2[i])*c - m2*s*((l1*c*(t1dot[i]**2))+l2*(t2dot[i]**2)) - (m1+m2)*g*np.sin(t1[i]))/(l1*(m1+m2*s*s))
        t2_acc[i] = (l1*s*(t1dot[i]**2) - g*np.sin(t2[i])-l1*c*t1_acc[i])/l2
    return t1_acc,t2_acc

def get_energy(res, cfg):
    
    '''Input - Position and Velocities at all time instances, length and mass and g
       Returns total energy at every time instance'''
    l1 = cfg.l1
    l2 = cfg.l2
    m1 = cfg.m1
    m2 = cfg.m2
    g = cfg.g
    
    t1 = res.yvec[:,0]
    t2 = res.yvec[:,1]
    t1dot,t2dot = get_velocities(res, cfg)
    
    T,V = [0]*len(t1),[0]*len(t1)
    for i in range(len(t1)):
        c = np.cos(t1[i] - t2[i])
        s = np.sin(t1[i] - t2[i])
        T[i] = (1/2)*m1*(l1**2)*(t1dot[i]**2) + (1/2)*m2*(((l1**2)*(t1dot[i]**2))+((l2**2)*(t2dot[i]**2)) + 2*l1*l2*t1dot[i]*t2dot[i]*c)
        V[i] = -m1*g*l1*np.cos(t1[i]) - m2*g*l1*np.cos(t1[i]) - m2*g*l2*np.cos(t2[i])
    return np.add(T,V)
    
def animate(res, cfg):
    filename = cfg.animTitle
    keepTrace = cfg.retainTrace
    
    logger.info('Animating '+filename+'.mp4')
    
    fig, ax = plt.subplots()
    maxsize=cfg.l1+cfg.l2
    ax.set(xlim=(-maxsize, maxsize), ylim=(-maxsize, maxsize))
    ax.set_aspect('equal')
    
    
    X1,Y1,X2,Y2 = get_positions(res, cfg)
    
    sizeCof = cfg.m2/cfg.m1;
    
    line = ax.plot([0, X1[0], X2[0]], [0, Y1[0], Y2[0]])[0]
    trace1 = ax.plot([0],[0], color='k', lw=1)[0]
    trace2 = ax.plot([0],[0], color='k', lw=1)[0]
    mass1 = ax.plot([X1[0]], [Y1[0]], 'o', color='tab:blue', label='Mass 1', markersize=5)[0]
    mass2 = ax.plot([X1[0]], [Y1[0]], 'o', color='tab:orange', label='Mass 2', markersize=5*sizeCof)[0]
    text = ax.text(-3.5, 3.5, '', fontsize=12)
        
    ax.set_xlabel('X, m')
    ax.set_ylabel('Y, m')
    ax.legend()
    
    def animate(i):
        line.set_xdata([0, X1[i], X2[i]])
        line.set_ydata([0, Y1[i], Y2[i]])
        mass1.set_xdata([X1[i]])
        mass1.set_ydata([Y1[i]])
        mass2.set_xdata([X2[i]])
        mass2.set_ydata([Y2[i]])
        if(keepTrace):
            iters = i
        else:
            iters = np.minimum(20,i)
        trace1X=np.zeros(iters)
        trace1Y=np.zeros(iters)
        trace2X=np.zeros(iters)
        trace2Y=np.zeros(iters)
        for k in range(0,iters):
            trace1X[k] = X1[i-k]
            trace1Y[k] = Y1[i-k]
            trace2X[k] = X2[i-k]
            trace2Y[k] = Y2[i-k]
        trace1.set_xdata(trace1X)
        trace1.set_ydata(trace1Y)
        trace2.set_xdata(trace2X)
        trace2.set_ydata(trace2Y)
        text.set_text(str(np.round(res.tvec[i],2))+" s")
    
    plt.rcParams['animation.ffmpeg_path'] = cfg.ffmpeg
    dt = cfg.dt
    fps = 1/dt
    anim = FuncAnimation(fig, animate, interval=1000/fps, frames=(res.tvec.size), repeat=True)
    # HTML(anim.to_html5_video()) #to play the animation in Jupyter Notebook
    if(cfg.animate):        
        plt.show()    
    if(cfg.mp4save):        
        anim.save(filename+'.mp4', fps=1/dt)
    