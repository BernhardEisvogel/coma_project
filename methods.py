#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:33:03 2020

@author: be
"""


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
    #konvertieren der Parameter zu floats
    for j in (0,len(data)-1):
        data[j]=float(data[j])
        
    #Bezeichnung der Parameter
    gamma,beta0,T,s0,i0=data[0],data[1],data[2],data[3],data[4]
    file.close()
    return gamma,beta0,T,s0,i0

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
    gamma, beta0,my,T, s0,i0 =data[0],data[1],data[2],data[3],data[4],data[5]
    file.close()
    return gamma, beta0, my, T, s0,i0