#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:35:43 2020

@author: Louisa Weber, Bernhard Eisvogel
"""
import methods
import numpy as np ## nump schon in Folgen importiert
import matplotlib.pyplot as plot

def Berechnung(x,y):
    return x,y-0.1

def Phasenportrait():
    xstartwerte= []
    ystartwerte = []
    for s in np.arange(0,1.05,0.1):
        xstartwerte.append(1-s)
        ystartwerte.append(s)     
    plot.plot(xstartwerte,ystartwerte ,color = 'b')
    for i in range(len(xstartwerte)):
        x= []
        y= []
        xw=xstartwerte[i]
        yw=ystartwerte[i]
        while yw>-0.000001:
            x.append(xw)
            y.append(yw)
            xw,yw = Berechnung(xw,yw)
        plot.plot(x,y,'b--')  
    plot.ylim(0, 1)
    plot.xlim(0, 1)
    plot.xlabel("Anteil Anfällige s")
    plot.ylabel("Anteil Anfällige I")
    plot.show
    return True

Phasenportrait()