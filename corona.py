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

def EpVerlauf():
    gamma,beta0,T,s0,i0 = methods.SirLesen()
    t=[]
    s=[]
    i=[]
    for r in np.arange(0,T+0.1,0.1):
        t.append(r)
        s.append(methods.epidlös(gamma,beta0, r,s0,i0,0)[0])
        i.append(methods.epidlös(gamma,beta0, r,s0,i0,0)[1])
    plot.plot(t,s)
    plot.plot(t,i)
    plot.title("Zeitlicher Verlauf der Anfälligen und Infizierten")
    plot.xlabel("Zeit t")
    plot.ylabel("Anteil der Anfälligen und Infizierten")
    plot.legend(["Anfällige s(t)", "Infizierte i(t)"])
    plot.show()

def Phasenportrait():
    gamma,beta0,T,s0,i0 = methods.SirLesen()
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
            xw,yw, a = methods.epidlös(gamma, beta0, 3,xw,yw,0)
            # besitzt keine weitere Bedeutung verhindert aber, 
            # dass eine zweite FUnktion geschrieben werden muss
        
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

gamma, beta0, my, T, s0,i0 = methods.SirDynLesen()
#my = 0.00003
#print("Infizierte in D nach", T, " Tagen: ", methods.endlös(my,gamma, beta0, T,s0,i0,0) )
#print("Infizierte in D nach", T, " Tagen: ", methods.epidlös(gamma, beta0, T,s0,i0,0) )
print(EpVerlauf())
#print(Phasenportrait())

