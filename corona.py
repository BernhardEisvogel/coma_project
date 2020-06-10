#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:35:43 2020

@author: Louisa Weber, Bernhard Eisvogel
"""
import methods
import numpy as np ## nump schon in Folgen importiert
import matplotlib.pyplot as plot
from matplotlib.lines import Line2D 

def Berechnung(x,y):
    return x,y-0.1

#def Phasenportrait():
xstartwerte= []
ystartwerte = []
for s in np.arange(0.0,1.05,0.1):
    xstartwerte.append(1-s) # Anfällige
    ystartwerte.append(s)   # Infektiös 
plot.plot(xstartwerte,ystartwerte ,color = 'b')
for i in range(len(xstartwerte)):
    x= np.array([])
    y= np.array([])
    xw=xstartwerte[i] 
    yw=ystartwerte[i]
    while yw>0.000001:
        x = np.append(x,[xw])
        y = np.append(y,[yw])
        xw,yw = methods.epidlös(3,xw,yw)
    
    plot.plot(x,y,'b--')
    if(len(x) != 0 and len(y) != 0):
        line = plot.plot(x,y,'b--')[0]
        methods.AddArrow(line, position = 0.05, color='b')
    
plot.ylim(0, 1)
plot.xlim(0, 1)
plot.title("Phasenportrait des SIR Modells")
plot.xlabel("Anteil Anfällige s")
plot.ylabel("Anteil Infektiöse i")
plot.show

#Phasenportrait()

