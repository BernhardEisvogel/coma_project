#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 13:56:28 2020

@author: be
"""

import math
# Dieses Programm enthält die Eulerapproximierung

# y0  = Anfangsbedinugng
# T   = Endzeit
# n   = Anzahl der Schritte
# fun = Die erste Ableitung der Funktion

# Hinweise: Es gibt n+1 Punkte, n muss später groß gewählt werden !

# Anfangszeit immer T=0 !!

def ForwardEuler(fun, y0, T,n):
    Schritt = y0           #Initialisieren
    for i in range(n+1): 
        Schritt = Schritt +  T/n*fun(T/n * i, Schritt)
    return Schritt

def ForwardKutta2_mod(fun, y0, T, n):
    # Diese Methode ist auch bekannt als modifizierte Euler'sche Methode
    Schritt = y0           #Initialisieren
    h = T/n
    for i in range(n+1): 
        Schritt = Schritt +  h* fun(h* i + 0.5 * h , Schritt + 
                                     0.5 * h * fun (h* i, Schritt))
    return Schritt

def ForwardKutta2_imp(fun, y0, T, n):
    # Diese Methode ist auch bekannt als verbesserte Euler'sche Methode
    Schritt = y0           #Initialisieren
    h = T/n
    for i in range(n+1): 
                Schritt = Schritt +  h * (0.5*fun(h*i,Schritt) + 
                                        0.5* fun (i*h + 0.5*h, Schritt 
                                                + h))
    return Schritt
    
def ForwardKutta3(fun, y0, T, n):
    # Runge Kutta Methode dritter Ordnung
    Schritt = y0           #Initialisieren
    h = T/n
    for i in range(n+1): 
        g1 = fun(h * i, Schritt)
        g2 =  fun(h*(i+0.5),Schritt + h/2 * g1)
        g3 =  fun(h*(i+0.5),Schritt + h/2 * g2)
        g4 =  fun(h*i,Schritt + h * g3)
        Schritt = Schritt +  T/n*(1/6 * g1+ 1/3 * g2 + 1/3 * g3 + 1/6 *g4)
    return Schritt




def TestF(x,y): return float(pow(math.e,x))

print("modifizierte Euker:", ForwardKutta2_mod(TestF, 0,6,1000))
print("verbesserte Euker:", ForwardKutta2_imp(TestF, 0,6,1000))
print("Euler:", ForwardEuler(TestF, 0,6,1000))
print("Kutta3:", ForwardKutta3(TestF, 0,6,1000))

print(pow(math.e,6)-1)

