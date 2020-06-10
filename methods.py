#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:33:03 2020

@author: Louisa Weber, Bernhard Eisvogel
"""
import numpy as np
import math



def ForwardEuler(fun, y0, T,n):
    '''Forward Euler nimmt als Parameter: \n
        fun, Ableitung der Funktion abhänjgig von x und y \n
        y0,  Anfangswert \n
        T, Endzeitpunkt \n
        n, Anzahl der Iterationsschritte
         \n
        Als Approximationsverfahren wird das Euelerverfahren benutzt
        '''
    n = n+1 # Es sollen n + 1 Punkte diskretisiert werden  
    Schritt = y0           
    for i in range(n): 
        Schritt = Schritt +  T/n*fun(T/n * i, Schritt)
    return Schritt

def ForwardKutta4(fun, y0, T, n):
    '''ForwardKutta4 nimmt als Parameter: \n
        fun, Ableitung der Funktion abhänjgig von x und y \n
        y0,  Anfangswert \n
        T, Endzeitpunkt \n
        n, Anzahl der Iterationsschritte \n
        
        Als Approximationsverfahren wird ein Vierstufiges Runge Kutta Verfahren
        angewendet
        '''    
    n = n+1                 # Es sollen n + 1 Punkte diskretisiert werden
    Schritt = y0           #Initialisieren
    h = T/n
    for i in range(n): 
        g1 = fun(h * i, Schritt)
        g2 =  fun(h*(i+1/2),Schritt + h * 1/2 * g1)
        g3 =  fun(h*(i+1/2),Schritt + h * 1/2 * g2)
        g4 =  fun(h * i + h,Schritt + h *       g3)
        Schritt = Schritt +  T/n*(1/6 * g1+ 1/3 * g2 + 1/3 * g3 + 1/6 *g4)
    return Schritt

def SirLesen():
    '''
    Dieser Funktion gibt das aus der Datei "sir.param" gelesen Tupel zurück
    -------
    gamma : FLOAT
        Erster Wert in [0]
    beta0 : FLOAT
        Zweiter Wert in [1]
    T : FLOAT
        Dritter Wert in [2]
    s0 : FLOAT
        Vierter Wert in [3]
    i0 : FLOAT
        Fünfter Wert in [4]
        '''
        
    file=open("sir.param","r")
    data=[]
    for i in file:
        #einzelne Parameter werden getrennt und in Liste gespeichert
        data+=i.split() 
    gamma=float(data[0])
    beta0=float(data[1])
    T=float(data[2])
    s0=float(data[3])
    i0=float(data[4])
    file.close()
    return gamma,beta0,T,s0,i0

gamma,beta0,T,s0,i0 = SirLesen()

def SirDynLesen():
    '''
    Dieser Funktion gibt das aus der Datei "sirdyn.param" gelesen Tupel zurück
    -------
    gamma : FLOAT
        Erster Wert in [0]
    beta0 : FLOAT
        Zweiter Wert in [1]
    my : FLOAT
        Dritter Wert in [2]
    T : FLOAT
        Vierter Wert in [3]
    s0 : FLOAT
        Fünfter Wert in [4]
    i0 : FLOAT
        Sechster Wert in [5]

    '''
    file=open("sirdyn.param","r")
    data=[]
    for i in file:
        #einzelne Parameter werden getrennt und in Liste gespeichert
        data+=i.split() 
    #konvertieren der Parameter zu floats
    for j in (0,len(data)-1):
        data[j]=float(data[j])
        
    #Bezeichnung der Parameter
    gamma, beta0, my, T, s0,i0 =data[0],data[1],data[2],data[3],data[4],data[5]
    file.close()
    return gamma, beta0, my, T, s0,i0

def ableitung(t,y):
    """
    berechnet einen Vektor der Änderunsraten von s und i

    Parameters
    ----------
    t : float, Zeitpunkt an dem die Ableitung berechnet wird
    y : array, y[1] Anzahl Anfälliger s(t)\n
    y[2] Anzahl Infizierter i(t)

    Returns
    -------
    Array, enthält Änderungsrate Anfälliger s und \n
    Änderungsrate Infizierter i

    """
    dS=(-1)*beta0*y[0]*y[1] #Ableitungsfunktion von s(t)
    dI=beta0*y[0]*y[1]-gamma*y[1] #Ableitungsfunktiono von i(t)
    return np.array([dS,dI])

def epidlös(t,s_0,i_0):
    """
    Lösung des epidemischen Modells unter Benutzung von Euler \n
    Funktion der Anfälligen s(t) \n
    und Infizierten i(t) 

    Parameters
    ----------
    t : float, Zeitpunkt  
    s0 : float, Anfangswerte Anfällige 
    i0 : float, Anfangswert Infizierte

    Returns
    -------
    array, enthält normierte Anzahl Anfällige s(t) \n
    und Infizierte i(t)

    """
    y_0=np.array([s_0,i_0])
    if t==0:
        return y_0
    else:
        return ForwardEuler(ableitung,y_0,t,1000)
    #n gross wählen!!!


