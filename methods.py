#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:33:03 2020

@author: Louisa Weber, Bernhard Eisvogel
"""
import numpy as np

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
    gamma=float(data[0])
    beta0=float(data[1])
    my=float(data[2])
    T=float(data[3])
    s0=float(data[4])
    i0=float(data[5])
    file.close()
    return gamma, beta0, my, T, s0,i0


def ForwardEuler(fun, y0, T,n):
    '''Forward Euler nimmt als Parameter: \n
        fun, Ableitung der Funktion abhängig von x und y \n
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
    Schritt = y0            #Initialisieren
    h = T/n
    
    for i in range(n): 
        g1 = fun(h * i, Schritt)
        g2 =  fun(h*(i+1/2),Schritt + h * 1/2 * g1)
        g3 =  fun(h*(i+1/2),Schritt + h * 1/2 * g2)
        g4 =  fun(h * i + h,Schritt + h *       g3)
        Schritt = Schritt +  T/n*(1/6 * g1+ 1/3 * g2 + 1/3 * g3 + 1/6 *g4)
    return Schritt




def epidlös(gamma, beta0, t,s_0,i_0,r_0):
    '''
    Lösung des epidemischen Modells unter Benutzung von Euler \n
    Funktion der Anfälligen s(t) \n
    und Infizierten i(t) 

    Parameters
    ----------
    gamma : float, Inverse der Durchschnittlichen Zeit, \n
    in der ein erkrankter ansteckend ist
    beta0 : float, Kontaktrate
    t : float, Zeitpunkt  
    s0 : float, Anfangswerte Anfällige 
    i0 : float, Anfangswert Infizierte

    Returns
    -------
    array, enthält normierte Anzahl Anfällige s(t) \n
    und Infizierte i(t)

    '''
    y_0=np.array([s_0,i_0,r_0])
    if t==0:
        return y_0
    else:
        return ForwardEuler(lambda t,y: 
                            np.array([(-1)*beta0*y[0]*y[1],
                                      beta0*y[0]*y[1]-gamma*y[1], gamma * y[1]]),
                                     y_0,
                                     t,
                                     2000)
    #n gross wählen!!!
            
           
def endlös(my, gamma, beta0, t,s_0,i_0, r_0):
    '''
    Lösung des epidemischen Modells unter Benutzung von Euler \n
    Funktion der Anfälligen s(t) \n
    und Infizierten i(t) 

    Parameters
    ----------
    N: Anzahl der Menschen
    gamma : float, Inverse der Durchschnittlichen Zeit, \n
    in der ein erkrankter ansteckend ist
    beta0 : float, Kontaktrate
    t : float, Zeitpunkt  
    s0 : float, Anfangswerte Anfällige 
    i0 : float, Anfangswert Infizierte

    Returns
    -------
    array, enthält Anzahl Anfällige s(t) \n
    und Infizierte i(t)

    '''
    y_0 = np.array([s_0,i_0,r_0])
    if t==0:
        return y_0
    else:
        return ForwardEuler(lambda t,y: 
                            np.array([my - my* y[0]-beta0*y[0]*y[1],
                                      (-1) * my * i_0 + beta0*y[0]*y[1]-gamma * y[1],
                                      (-1) * my*y[2]+ gamma * y[1]]),
                                     y_0,
                                     t,
                                     2000)
    #n gross wählen!!!


"""
Created on Thu Jun  4 06:30:03 2020

@author: Lena
"""


def AddArrow(line, position=None, direction='right', size=15, color=None):
    '''
    Add an arrow to a line.

    line:       Line2D object
    position:   y-position of the arrow. If None, mean of ydata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    '''
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = ydata.mean()
        
    # find closest index
    start_ind = np.argmin(np.abs(ydata - position))
    
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )

