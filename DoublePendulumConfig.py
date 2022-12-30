# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 18:02:15 2020

@author: Group 5
"""

import xml.etree.ElementTree as xmlET
import logging
import numpy as np

logger = logging.getLogger()

class DoublePendulumConfig:
    solver = 1
    endtime = 10
    dt = 0.01
    
    positionPlot = False
    energyPlot = False
    velocityPlot = False
    accelerationPlot = False
    animate = False
    mp4save = False
    ffmpeg = ""
    animTitle = ""
    retainTrace = True
    
    m1 = 1
    m2 = 1
    l1 = 1
    l2 = 1
    
    t1 = 0
    t2 = 0
    p1 = 0
    p2 = 0

def xmlprint(xmlele, level):
    spacer = "";
    for k in range(0,level):
        spacer = spacer.join("  ")
    if(len(xmlele.getchildren())==0):
        content = spacer + xmlele.tag + " " + xmlele.text
        logger.debug(content)
    else:
        content = spacer + xmlele.tag
        logger.debug(content)
    
    for xmlchild in xmlele:
        xmlprint(xmlchild, level+1)

def readConfig(cfg):
    xmltree = xmlET.parse(cfg)
    xmlroot = xmltree.getroot()
    logger.debug("root name: "+xmlroot.tag)
    # xmlprint(xmlroot, 0)
    
    config = DoublePendulumConfig()
    readProblemElement(xmlroot.find('problem'), config)
    readSimulationElement(xmlroot.find('sim'), config)
    readPostProcessingElement(xmlroot.find('postprocessing'), config)
    return config
    
def readProblemElement(element, cfg):    
    if(element is not None):
        cfg.g = getFloatIfNotEmpty(element, 'g', 9.81)
        for c in element.getchildren():
            if(c.tag == 'body'):
                readBodyElement(c, cfg)
        logger.info('Problem definition read')
    else:
        logger.warn("Problem element not found in config, using default values")

def readBodyElement(element, cfg):
    eleid = getIntParameterIfExists(element, 'id', 0)
    if(eleid == 1):
        cfg.m1 = getFloatIfNotEmpty(element, 'm', 1)
        cfg.l1 = getFloatIfNotEmpty(element, 'l', 1)
        cfg.t1 = (getFloatIfNotEmpty(element,'theta', 45))*np.pi/180
        cfg.p1 = getFloatIfNotEmpty(element, 'p', 0)
    elif(eleid == 2):
        cfg.m2 = getFloatIfNotEmpty(element, 'm', 1)
        cfg.l2 = getFloatIfNotEmpty(element, 'l', 1)
        cfg.t2 = (getFloatIfNotEmpty(element,'theta', 90))*np.pi/180
        cfg.p2 = getFloatIfNotEmpty(element, 'p', 0)
        
def readSimulationElement(element, cfg):    
    if(element is not None):  
        cfg.solver = getIntIfNotEmpty(element, 'solver', 3)
        cfg.dt = getFloatIfNotEmpty(element, 'dt', 3)
        cfg.endtime = getFloatIfNotEmpty(element, 'endtime', 3)
        logger.info('Simulation definition read')
    else:
        logger.warn("Simulation element not found in config, using default values")

def readPostProcessingElement(element, cfg):    
    if(element is not None):    
        cfg.positionPlot = getBoolParameterIfExists(element.find('positionplot'), 'draw', False)
        cfg.velocityPlot = getBoolParameterIfExists(element.find('velocityplot'), 'draw', False)
        cfg.accelerationPlot = getBoolParameterIfExists(element.find('accelerationplot'), 'draw', False)
        cfg.energyPlot = getBoolParameterIfExists(element.find('energyplot'), 'draw', False)
        cfg.animate = getBoolParameterIfExists(element.find('animation'), 'draw', False)
        cfg.mp4save = getBoolParameterIfExists(element.find('animation'), 'mp4save', False)
        cfg.retainTrace = getBoolParameterIfExists(element.find('animation'), 'retainTrace', False)
        cfg.ffmpeg = getStringParameterIfExists(element.find('animation'), 'ffmpeg', '')
        cfg.animTitle = getStringParameterIfExists(element.find('animation'), 'title', 'animation')
        logger.info('Post processing definition read')
    else:
        logger.warn("PostProcessing element not found in config, using default values")        

def getIntParameterIfExists(element, key, default):
    if(element is not None):
        try:
            value = element.attrib[key]
            intval = int(value)
            return intval
        except:
            logger.warn("Element "+element.tag+" does not have an int attribute "+key)
            return default;
    else:
        return default

def getStringParameterIfExists(element, key, default):
    if(element is not None):
        try:
            value = element.attrib[key]
            string = str(value)
            return string
        except:
            logger.warn("Element "+element.tag+" does not have a string attribute "+key)
            return default;
    else:
        return default
    
def getBoolParameterIfExists(element, key, default):
    if(element is not None):
        try:
            value = element.attrib[key]
            intval = int(value)
            boolean = bool(intval)
            return boolean
        except:
            logger.warn("Element "+element.tag+" does not have a boolean attribute "+key)
            return default;
    else:
        return default

def getFloatIfNotEmpty(parentXml, tag, default):
    child = parentXml.find(tag)
    if(child is not None):
        try :  
            fl = float(child.text) 
            return fl
        except : 
            logger.warn("Element "+tag+" of parent "+parentXml.tag+" does not have a float")
            return default
    else:
        return default
    
def getIntIfNotEmpty(parentXml, tag, default):
    child = parentXml.find(tag)
    if(child is not None):
        try :  
            fl = int(child.text) 
            return fl
        except : 
            logger.warn("Element "+tag+" of parent "+parentXml.tag+" does not have an int")
            return default
    else:
        return default
    



        