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




def epidlös(gamma, beta0, t,s_0,i_0,r_0, schritte = 2000, kutta=False):
    '''
    Lösung des epidemischen Modells unter Benutzung von Euler/RungeKutta \n
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
    schritte : int, Dieser Parameter ist optional (standartwert 2000) 
    und gibt die Anzahl der Schritte an, die bis zum Zeitpunkt t durchgeführt 
    werden sollen
    kutta : bool, dieser Parameter ist optional (standartwert false) und gibt an, 
    ob an Stelle des Eulerverfahrens das Rungeuuttaverfahren benutzt werden soll.

    Returns
    -------
    array, enthält normierte Anzahl Anfällige s(t), \n
    Infizierte i(t) und genesene r(t)

    '''
    y_0=np.array([s_0,i_0,r_0])
    if t==0:
        return y_0
    else:
        if(kutta==False):
            return ForwardEuler(lambda t,y: 
                                np.array([(-1)*beta0*y[0]*y[1],
                                          beta0*y[0]*y[1]-gamma*y[1], gamma * y[1]]),
                                         y_0,
                                         t,
                                         schritte)
        else: 
            return ForwardKutta4(lambda t,y: 
                                np.array([(-1)*beta0*y[0]*y[1],
                                          beta0*y[0]*y[1]-gamma*y[1], gamma * y[1]]),
                                         y_0,
                                         t,
                                         schritte)
    #n gross wählen!!!
            
           
def endlös(my, gamma, beta0, t, s_0, i_0, r_0, schritte = 2000):
    '''
    Lösung des endemischen Modells unter Benutzung von Euler/Rungekutta \n
    Funktion der Anfälligen s(t) \n
    und Infizierten i(t) 

       Parameters
    ----------
    my : float, Geburtenrate
    gamma : float, Inverse der Durchschnittlichen Zeit, \n
    in der ein erkrankter ansteckend ist
    beta0 : float, Kontaktrate
    t : float, Zeitpunkt  
    s0 : float, Anfangswerte Anfällige 
    i0 : float, Anfangswert Infizierte
    schritte : int, Dieser Parameter ist optional (standartwert 2000) 
    und gibt die Anzahl der Schritte an, die bis zum Zeitpunkt t durchgeführt 
    werden sollen

    Returns
    -------
    array, enthält normierte Anzahl Anfällige s(t), \n
    Infizierte i(t) und genesene r(t)
'''
    y_0 = np.array([s_0,i_0,r_0])
    if t==0:
        return y_0
    else:
        return ForwardEuler(lambda t,y: 
                            np.array([my - my* y[0]-beta0*y[0]*y[1],
                                      (-1) * my * y[1] + beta0*y[0]*y[1]-gamma * y[1],
                                      (-1) * my*y[2]+ gamma * y[1]]),
                                     y_0,
                                     t,
                                     schritte)
    #n gross wählen!!!
         
def TabelleLesen():
    """
    Rückgabe der aus der Datei "Tabellendaten.txt" eingelesenen Daten 

    Returns
    -------
    list, enthält Infiziertendaten in Deutschland vom 1.3. bis 27.5.

    """
    file=open("tabellendaten.txt","r")
    inf=[]
    for i in file:
        inf+=i.split()
    for i in range(0,len(inf),1):
        inf[i]=int(inf[i]) 
    return inf
    
def fehler(t,beta,y):
    """
    Berechnung des Fehler zwischen der errechneten \n
    Infiziertenzahl des endemischen Modells und den Tabellendaten \n
    zum Zeitpunkt t in Abhängigkeit von der gewählten Kontaktrate und \n
    den Anfangswerten y
    Parameters
    ----------
    t : integer, Zeitpunkt 
    beta : float, Kontaktrate
    y: list, enthält s0,i0,r0 zur Berechnung der Infiziertenzahl mit endlös()
    Returns
    -------
    fehler : relativer Fehler zwischen Berechnung und tatsächlicher Infiziertenzahl
    (ohne Absolutbetrag)
    """
    inf=TabelleLesen()
    N=83*(10**6)
    my=1/27375
    gamma=1/6.5
    fehler=(endlös(my,gamma,beta,t,y[0],y[1],y[2])[1]*N-inf[t])/inf[t]
    return fehler

def betaapprox(t,beta0,a,y):
    """
    Approximation der Kontaktrate für einen Zeitpunkt t \n
    sodass der relative Fehler zwischen berechneter und aktueller \n 
    Infiziertenzahl kleiner 0.05 ist, bei einem Anfangswert von beta0 \n
    und den Anfangswerten y
    Parameters
    ----------
    t : integer, Zeitpunkt
    beta0 : float, Anfangswert der Approximation der Kontaktrate
    a : integer, Wert -1 oder 1, Richtung der Approximation von beta0 ausgehend
    y: list, enthält s0,i0,r0 
    Returns
    -------
    beta : float, approximierte Kontaktrate
    """
    beta=beta0
    eps=abs(fehler(t,beta,y))#relativer Fehler im Betrag beim Anfangswert beta0
    while (eps>=0.05 and beta<1):
        beta=round(beta+(a*0.001),4)
        eps=abs(fehler(t,beta,y))
    #Annäherung von beta, solange bis Fehler kleiner 0.05
    return beta

def Kontaktrate(t):
    """
    stückweise konstante Funktion der Kontaktrate bis zum Zeitpunkt t \n
    unter Verwendung von betaapprox, fehler
    Parameters
    ----------
    t : integer, Zeitpunkt
    Returns
    -------
    list, Kontaktrate zu den Zeitpunkten 0,1,...,t
    """
    N=83*(10**6)
    my=1/27375
    gamma=1/6.5
    time=0
    beta=betaapprox(1,0,1,[(N-130)/N,130/N,0])
    #approximierte Kontaktrate bei t=1, beta0=0, da bei t=0 keine sinnvolle \n
    #Approximation stattfindet, da der Fehler zum Anfangswert i0 hier immer 0
    eps=0
    y=[(N-130)/N,130/N,0]
    kontakt=np.array([])
    while time<t+1:
        while (eps<0.05 and time<t+1):
            kontakt=np.append(kontakt,[beta])
            y=endlös(my,gamma,beta,time,y[0],y[1],y[2])
            time+=1
            if time<t+1:
                eps=abs(fehler(time,beta,y))
            #Kontaktrate konstant lassen, bis der Fehler größer 0.05 \n
        if time<t+1:
            #Neue Approximation der Kontaktrate falls Fehler größer 0.05
            a=(-1)*fehler(time,beta,y)/abs(fehler(time,beta,y))
            beta=betaapprox(time,beta,a,y)
            eps=abs(fehler(time,beta,y))
    return kontakt[0:(t+1)]
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

def maximalWert(gamma, beta0, t, s_0, i_0, r_0, p, schritte = 2000):
    '''
    Lösung des epidemischen Modells unter Berücksichtung der Impfungsrate p.
    Gibt außerdem max(i) zurück

    Parameters
    ----------
    gamma : float, Inverse der Durchschnittlichen Zeit, \n
    in der ein erkrankter ansteckend ist
    beta0 : float, Kontaktrate
    t : float, Zeitpunkt  
    s0 : float, Anfangswerte Anfällige 
    i0 : float, Anfangswert Infizierte
    schritte : int, Dieser Parameter ist optional (standartwert 2000) 
    und gibt die Anzahl der Schritte an, die bis zum Zeitpunkt t durchgeführt 
    werden sollen
    kutta : bool, dieser Parameter ist optional (standartwert false) und gibt an, 
    ob an Stelle des Eulerverfahrens das Rungeuuttaverfahren benutzt werden soll.

    Returns
    -------
    array, enthält normierte Anzahl Anfällige s(t), \n
    Infizierte i(t) und genesene r(t)\
    
    float, maximum der Infizierten

    '''

    y_0=np.array([s_0,i_0,r_0])
    if t==0:
        return y_0
    else:
        return ForwardEulerMitMax(lambda t,y: 
                            np.array([(-1)*beta0*y[0]*y[1] - p,
                                      beta0*y[0]*y[1]-gamma*y[1], 
                                      gamma * y[1] +p]),
                                     y_0,
                                     t,
                                     schritte)

    #n gross wählen!!!
def ForwardEulerMitMax(fun, y0, T,n):
    '''ForwardEulerMitMax nimmt als Parameter: \n
        fun, Ableitung der Funktion abhängig von x und y \n
        y0,  Anfangswert \n
        T, Endzeitpunkt \n
        n, Anzahl der Iterationsschritte
         \n
        Als Approximationsverfahren wird das Euelerverfahren benutzt
        
        Zusätzlich zur Lösung gibt diese Funktion auch das Maximum des zweiten
        Wertes des Eingabevektors zurück.
        '''
    n = n+1 # Es sollen n + 1 Punkte diskretisiert werden  
    Schritt = y0         
    maximalInfiziert = 0
    for i in range(n): 
        Schritt = Schritt +  T/n * fun(T/n * i, Schritt)
        if(Schritt[1] > maximalInfiziert):
            maximalInfiziert = Schritt[1]
        Schritt[1] = identiodernull(Schritt[1])
        Schritt[0] = identiodernull(Schritt[0])
    return Schritt, maximalInfiziert

def identiodernull(r):
    if(r<0): return 0
    else: 
        return r
