# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:59:02 2020

@author: Group 5
"""

import argparse
import logging
import numpy as np

from DoublePendulumConfig import readConfig
from DoublePendulum import DoublePendulum
from DoublePendulum import plot_positions
from DoublePendulum import plot_velocities
from DoublePendulum import plot_accelerations
from DoublePendulum import plot_energy
from DoublePendulum import animate

from ImplicitEulerNewton import ImplicitEulerNewton
from DefaultSchemes import ExplicitEuler
from DefaultSchemes import ImplicitEuler
from DefaultSchemes import Midpoint
from DefaultSchemes import RK4
import Simulationresults

logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s", level=logging.INFO)

logger = logging.getLogger()

argparser = argparse.ArgumentParser()
argparser.add_argument("config", type=str)
argparser.add_argument("-v", "--verbose", action="count")
args = argparser.parse_args()

if args.verbose is not None:
    if args.verbose>1:
        logger.setLevel(logging.DEBUG)
        logging.info("LOG level set to DEBUG")
    else:
        logger.setLevel(logging.INFO)
        logging.info("LOG level set to INFO")
        
cfg = readConfig(args.config)

def RHS(t, x):
    return DoublePendulum(t,x,cfg)

initialCond = np.array([cfg.t1, cfg.t2, cfg.p1, cfg.p2])
logger.debug('Initial conditions: '+' '.join(str(initialCond)))

solver_type = cfg.solver
end = cfg.endtime
dt = cfg.dt

if(solver_type==1):
    logging.info("Solving with Explicit Euler")
    res = ExplicitEuler(RHS, 0, end, initialCond, dt)
elif(solver_type==2):
    logging.info("Solving with Implicit Euler /w Newton-Raphson scheme")
    res = ImplicitEulerNewton(dt, end, cfg, initialCond)
elif(solver_type==3):
    logging.info("Solving with Explicit Classical 4th order Runge Kutta method")
    res = RK4(RHS, 0, end, initialCond, dt)
elif(solver_type==4):
    logging.info("Solving with Implicit Euler /w fixed point iteration")
    res = ImplicitEuler(RHS, 0, end, initialCond, dt)
elif(solver_type==5):
    logging.info("Solving with Midpoint rule /w fixed point iteration")
    res = Midpoint(RHS, 0, end, initialCond, dt)
else:
    logging.info("Solving with default solver RK4")
    res = RK4(RHS, 0, end, initialCond, dt)

if(cfg.positionPlot):
    logging.info("Plotting positions")
    plot_positions(res, cfg)
if(cfg.velocityPlot):
    logging.info("Plotting velocities")
    plot_velocities(res, cfg)
if(cfg.accelerationPlot):
    logging.info("Plotting accelerations")
    plot_accelerations(res, cfg)
if(cfg.energyPlot):
    logging.info("Plotting energy")
    plot_energy(res, cfg)
if(cfg.animate or cfg.mp4save):
    animate(res, cfg)
    
input("Press Enter to exit...") #if you run in Jupyter, press Enter to view the results