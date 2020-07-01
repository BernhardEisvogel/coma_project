#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:35:43 2020
@author: Louisa Weber, Bernhard Eisvogel
"""
import methods
import numpy as np
import matplotlib.pyplot as plot
import math
from scipy.stats import linregress 

def EpVerlauf(end = True):
    """
    Dieses Programm liest die Datei sirdyn.param und zeichnet für die daraus\n
    ausgelesenen Werte ein Diagramm mit dem Verlauf der Infizierten und An- \n
    fälligen Zahl
    
    Es wird standardmäßig das epidemologische Modell benutzt. Dies kann aber \n
    end=False als Eingabewert umgestellt werden.
    """
    gamma, beta0, my, Tgelesen, s0,i0 = methods.SirDynLesen()
    t=[0]
    s=[s0]
    i=[i0]
    rec=[0] # Dieser Wert wird zwar aktuell nicht benötigt, steht aber da für 
            # zukünftige Überlegungen
    if(end == False): my = 0
    
    for r in np.arange(0,Tgelesen+0.1,2):
        t.append(r)
        loesung = methods.endlös(my, gamma, beta0, 2, s[len(s)-1],i[len(i)-1],
                                 rec[len(rec)-1])
        s.append(loesung[0])
        i.append(loesung[1])
        rec.append(loesung[2])
    plot.plot(t,s)
    plot.plot(t,i)
    plot.title("Zeitlicher Verlauf der Anfälligen und Infizierten")
    plot.xlabel("Zeit t")
    plot.ylabel("Anteil der Anfälligen und Infizierten")
    plot.legend(["Anfällige s(t)", "Infizierte i(t)"])
    plot.show()
    plot.show()
def DatenSir():
    gamma,beta0,T,s0,i0 = methods.SirLesen()
    N=83*(10**6)
    t=np.array(["Zeitpunkt"])
    e=np.array(["Gesamtanzahl der Erkrankten"])
    r=np.array(["Gesamtanzahl der Genesenen"])
    for i in np.arange(0,T+0.1,1):
        t=np.append(t,[str(int(i))])
        e=np.append(e,[str(int((1-methods.epidlös(gamma,beta0,i,s0,i0,0)[0])*N))])
        r=np.append(r,[str(int(methods.epidlös(gamma,beta0,i,s0,i0,0)[2]*N))])
    daten=np.vstack((t,e,r))
    daten=daten.T
    np.savetxt("daten.sir.txt",daten,fmt="%20s")

def PhasenportraitEp():
    '''
    Dieses Programm zeichnet ein Phasenportrait des epidemologischen Epidemie-\n
    verlauf. Dazu liest es die Daten aus sirdyn.param aus

    Returns
    -------
    None.

    '''
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
        
        plot.plot(x,y,'b--')
        if(len(x) != 0 and len(y) != 0):
            line = plot.plot(x,y,'b--')[0]
            methods.AddArrow(line, position = 0.05, color='b')
        
    plot.ylim(0, 1)
    plot.xlim(0, 1)
    plot.title("Phasenportrait des epidemischen SIR Modells")
    plot.xlabel("Anteil Anfällige s")
    plot.ylabel("Anteil Infektiöse i")
    plot.show
    
def PhasenportraitEnd():
    '''
    Dieses Programm zeichnet ein Phasenportrait des endemischen Epidemie-\n
    verlauf mit den Daten aus sirdyn.param

    Returns
    -------
    None.

    '''
    gamma, beta0, my, T, s0,i0 = methods.SirDynLesen()
    #my=0.00003653
    xstartwerte=[]
    ystartwerte=[]
    for s in np.arange(0.0,1.05,0.1):
        xstartwerte.append(1-s)
        ystartwerte.append(s)
    plot.plot(xstartwerte,ystartwerte,color="b")
    for i in range(len(xstartwerte)):
        x=np.array([])
        y=np.array([])
        xw=xstartwerte[i]
        yw=ystartwerte[i]
        while yw>0.003:
            x=np.append(x,[xw])
            y=np.append(y,[yw])
            xw,yw,a=methods.endlös(my,gamma,beta0,2,xw,yw,0)
        plot.plot(x,y,"b--")
        if len(x)!=0:
            line=plot.plot(x,y,"b--")[0]
            methods.AddArrow(line, position=0.065,color="b")
    
    plot.ylim(0,1)
    plot.xlim(0,1)
    plot.title("Phasenportrait des endemischen SIR Modells")
    plot.xlabel("Anteil Anfällige s")
    plot.ylabel("Anteil Infektiöse i")
    plot.show()

def Verlaufaktuell():
    inf=methods.TabelleLesen()
    Kontakt=methods.Kontaktrate(len(inf)-1)
    t=[]
    I=[]
    N=83*(10**6)
    for i in range(0,len(inf),1):
        t.append(i)
        I.append(N*(methods.endlös(1/27375,1/6.5,Kontakt[i],i,(N-130)/N,
                                   130/N,0)[1]))

    plot.plot(t,inf)
    plot.plot(t,I)
    plot.title("Zeitlicher Verlauf der Infiziertenzahl in Deutschland")
    plot.xlabel("Datum")
    plot.ylabel("Anzahl der Infizierten")
    plot.xticks( [0, 30, 60, 87],
            [ "1.3.20", "31.03.20", "30.04.20", "27.05.20" ] )
    plot.legend(["Tabellendaten", "Errechnete Daten"])
    plot.show()
    
def KontaktrateVerlauf():
    x=[]
    y=[]
    k=methods.Kontaktrate(87)
    for i in range(0,88,1):
        x.append(i)
        y.append(k[i])
    plot.plot(x,y)
    plot.title("Zeitlicher Verlauf der Kontaktrate")
    plot.xlabel("Datum")
    plot.ylabel("Kontaktrate")
    plot.xticks( [0, 30, 60, 87],
            [ "1.3.20", "31.03.20", "30.04.20", "27.05.20" ] )
    plot.show()

def Prognose():
    """
    Prognose für die Infiziertenzahl am 1.7. ausgehend von \n
    der Infiziertenzahl am 24.6. unter der Annahme einer konstanten \n
    Kontaktrate ab 24.6.
    """
    N=83*(10**6)
    my=1/27375
    gamma=1/6.5
    beta=0
    datum=122 # Prognose für den 1.7. (122 Tage nach 1.3.)
    aktdat=115 # Prognose vom 24.6. (115 Tagen nach 1.3.)
    aktinf=6218 # aktuell Infizierte
    #Approximation der Kontaktrate vom 24.6. bis Fehler zwischen
    #berechneter und tatsächlicher Infiziertenzahl am 24.6. kleiner 0.05
    eps=abs((methods.endlös(my,gamma,beta,aktdat,(N-130)/N,130/N,0)[1]*N-aktinf)/aktinf)
    while (eps>=0.05 and beta<1):
        beta=round(beta+0.001,4)#evt kleiner wählen
        eps=abs((methods.endlös(my,gamma,beta,aktdat,(N-130)/N,
                                130/N,0)[1]*N-aktinf)/aktinf)
    #Progonse für 1.7. bei Kontaktrate vom 24.6.
    prognose=int(methods.endlös(my,gamma,beta,datum,(N-130)/N,130/N,0)[1]*N)
    return prognose


def fehlerBerechnenAbsolut(t):
    '''
    Diese Funktion zeichnet ein Schaubild zur Visualisierung der Ungenauigkeit \n
    des Eulerapproximationsverfahren. Nimmt als Variable den Zeitpunkt 0<=t<10\n
    an, für den der Fehler visualisiert werden soll
    Returns
    -------
    None.
    '''
    
    if(t<=0 or t>10):
        print("Bitte geben sie einen gültigen t-Wert ein!")
        exit()
    
    beta0 = 1
    gamma = 0
    s0   = -9
    i0    = 10
    r0 = 0
    x =[]
    i = []
    s = []
    echteWerte = [1-(math.exp(t))/(math.exp(t) -0.9),(math.exp(t))/(math.exp(t) -0.9),0]
    for n in 1000, 2000, 4000, 8000, 16000:
        approximiert = methods.epidlös(gamma, beta0, t, s0,i0,r0, schritte = n)
        x.append(n)
        i.append(abs(approximiert[1]-echteWerte[1]))
        s.append(abs(approximiert[0]-echteWerte[0]))
    plot.plot(x,i, 'b--',linewidth=2)
    plot.plot(x,s,'r-',linewidth=1)
    plot.title("Vergleich der Fehler der Approximierung")
    plot.xlabel("Anzahl der Schritte")
    plot.ylabel("Absoluter Fehler")
    plot.legend(["Fehler bei den Infizierte", "Fehler bei den Anfälligen"])
    plot.show()
        
def fehlerVergleichloglog(t):
    '''
    Diese Funktion zeigt die Ungenauigkeit des Euelrapproximationsverfahren auf
    einer doppelt Logarithmsichen Skala an und vergleicht ihn mit einer Funktion
    der Form f(n) = c/n zur Fehlervorhersage .Nimmt als Variable den Zeitpunkt 0<=t<10\n
    an, für den der Fehler visualisiert werden soll
    Returns
    -------
    None.
    '''
    if(t<=0 or t>10):
        print("Bitte geben sie einen gültigen t-Wert ein!")
        exit()
    
    
    beta0 = 1
    gamma = 0
    s0   = -9
    i0    = 10
    r0 = 0
    x =[]
    i = []
    s = []
    f= []
    echteWerte = [1-(math.exp(t))/(math.exp(t) -0.9),(math.exp(t))/(math.exp(t) -0.9),0]
    for n in 1000, 2000, 4000, 8000, 16000:
        approximiert = methods.epidlös(gamma, beta0, t, s0,i0,r0, schritte = n)
        x.append(n)
        i.append(abs(approximiert[1]-echteWerte[1]))
        s.append(abs(approximiert[0]-echteWerte[0]))
    plot.loglog(x,i, 'b--',linewidth=4)
    plot.loglog(x,s,'b-',linewidth=2)
    
    for j in range(0,len(i)):
        i[j] = math.log(i[j])
    xinlog = []
    for r in range(0,len(x)):
        xinlog.append(math.log(x[r]))
            
    print(linregress(xinlog,i))
    c = pow(math.e, linregress(xinlog,i).intercept) 
    for n in 1000, 2000, 4000, 8000, 16000:
        f.append(c/n)
    plot.loglog(x,f,'r-',linewidth=1)
    plot.title("Vergleich des Fehlers des Eulerverfahrens \n mit der"+ 
               "Fehlerfunktion f für t=" + str(t))
    plot.xlabel("Anzahl der Schritte")
    plot.ylabel("Absoluter Fehler")
    plot.legend(["Fehler bei den Infizierte", "Fehler bei den Anfälligen", 
                 "f-Funktion mit c= %.3f" % c])

def EpVerlaufMitImpfung(TBeginnImpfung = 6, TEndeImpfung = 80,p = 0.001):
    '''
    Diese Funktion zeigt die Ungenauigkeit des Euelrapproximationsverfahren auf
    einer doppelt Logarithmsichen Skala an und vergleicht ihn mit einer Funktion
    der Form f(n) = c/n zur Fehlervorhersage .Nimmt als Variable den Zeitpunkt 0<=t<10\n
    an, für den der Fehler visualisiert werden soll
    Returns
    -------
    None.
    '''
    
    TEndeModell = 130
    if(TEndeImpfung > TEndeModell): print("Sie haben ungültige Parameter für die Impf-\n" + 
                                          "dauer angegeben")
    beta0 = 0.472
    Bevoelkerung = 82000000
    maximaleInfiziertenanzahl = 350000
    gamma = 1/6.5
    i0=1/ 82000000
    s0=1-i0
    t=[0]
    s=[s0]
    i=[i0]
    rec=[0]
    for r in np.arange(0,TBeginnImpfung,1): #Vor der Impfung
        t.append(r)
        loesung = methods.endlös(0, gamma, beta0, 1, s[len(s)-1],i[len(i)-1],
                                 rec[len(rec)-1])
        s.append(loesung[0])
        i.append(loesung[1])
        rec.append(loesung[2])
    maximal = 0 # benötigt zur Betrachtung der Belastungsgrenze des Krankenhauses
    for r in np.arange(TBeginnImpfung,TEndeImpfung,1): #während der Impfperiode
        t.append(r)
        loesung, maximalneu = methods.maximalWert(gamma, beta0, 1, s[len(s)-1],
                                                  i[len(i)-1],rec[len(rec)-1],p)
        if(maximalneu > maximal): maximal = maximalneu
        s.append(loesung[0])
        i.append(loesung[1])
        rec.append(loesung[2])
    for r in np.arange(TEndeImpfung,TEndeModell,1): # nach der Impfperiode
        t.append(r)
        loesung, maximalneu = methods.maximalWert(gamma, beta0, 1, s[len(s)-1],
                                                  i[len(i)-1],rec[len(rec)-1],0)
        if(maximalneu > maximal): maximal = maximalneu
        s.append(loesung[0])
        i.append(loesung[1])
        rec.append(loesung[2])
        
    
    plot.plot(t,s,color='b')
    plot.plot(t,i,color='orange') 
    plot.legend(["Anfällige s(t)", "Infizierte i(t)"])
    plot.axvline(x=TBeginnImpfung-1, linewidth=0.5, color='r')
    plot.axhline(y=maximaleInfiziertenanzahl /Bevoelkerung, linewidth=0.5, color='r')
    plot.axvline(x=TEndeImpfung-1, linewidth=0.5, color='r')
    print("Die maximale Infiziertenanzahl beträgt: %.0f" %(maximal * Bevoelkerung))
    plot.title("Zeitlicher Verlauf der Anfälligen und Infizierten mit Impfung")
    plot.xlabel("Zeit t in Tagen")
    plot.ylabel("Anteil der Anfälligen und Infizierten")
    plot.show()
    plot.show()