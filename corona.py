#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:35:43 2020

@author: Louisa Weber, Bernhard Eisvogel
"""
import methods
import numpy as np ## nump schon in Folgen importiert
import matplotlib.pyplot as plot
#from matplotlib.lines import Line2D 

def EpVerlauf(end = 'False', alpha = 0, T = None):
    gamma, beta0, my, Tgelesen, s0,i0 = methods.SirDynLesen()
    t=[]
    s=[]
    i=[]
    if(T != None):
        Tgelesen = T
        
    if(end == 'False'):
        my = 0
        
    for r in np.arange(0,Tgelesen+0.1,2):
        t.append(r)
        s.append(methods.epidlösGesamt(my, gamma,beta0,alpha, r,s0,i0,0)[0])
        i.append(methods.epidlösGesamt(my, gamma,beta0,alpha, r,s0,i0,0)[1])
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
    plot.title("Phasenportrait des epidemischen SIR Modells")
    plot.xlabel("Anteil Anfällige s")
    plot.ylabel("Anteil Infektiöse i")
    plot.show
    
def PhasenportraitEnd():
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
    
def PhasenportraitFallBack(alpha):
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
            xw,yw,a=methods.epidlösGesamt(my,gamma,beta0,alpha,4,xw,yw,0)
        plot.plot(x,y,"b--")
        if len(x)!=0:
            line=plot.plot(x,y,"b--")[0]
            methods.AddArrow(line, position=0.065,color="b")
    
    plot.ylim(0,1)
    plot.xlim(0,1)
    plot.title("Phasenportrait des endemischen SIR Modells mit Fall Back")
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
        I.append(N*(methods.endlös(1/27375,1/6.5,Kontakt[i],i,(N-130)/N,130/N,0)[1]))

    plot.plot(t,inf)
    plot.plot(t,I)
    plot.title("Zeitlicher Verlauf der Infiziertenzahl in Deutschland")
    plot.xlabel("Datum")
    plot.ylabel("Anzahl der Infizierten")
    plot.xticks( [0, 30, 60, 87],
            [ "1.3.20", "31.03.20", "30.04.20", "27.05.20" ] )
    plot.legend(["Tabellendaten", "Errechnete Daten"])
    plot.show()


def Prognose():
    """
    Prognose für die Infiziertenzahl am 1.7. ausgehend von \n
    der Infiziertenzahl am 24.6. unter der Annahme einer konstanten \n
    Kontaktrate ab 24.6.

    Returns
    -------
    prognose : integer, prognostizierte Infiziertenzahl am 1.7.

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
        eps=abs((methods.endlös(my,gamma,beta,aktdat,(N-130)/N,130/N,0)[1]*N-aktinf)/aktinf)
    #Progonse für 1.7. bei Kontaktrate vom 24.6.
    prognose=int(methods.endlös(my,gamma,beta,datum,(N-130)/N,130/N,0)[1]*N)
    return prognose
#print(Prognose())
#print(Verlaufaktuell())

#PhasenportraitFallBack(0.25)
#EpVerlaufFallBack(0.25)
#EpVerlauf()
PhasenportraitFallBack(0.25)

# print("listo")
# gamma, beta0, my, T, s0,i0 = methods.SirDynLesen()
# T = 200
# print("End, Infizierte in D nach", T, " Tagen: ", methods.endlös(my,gamma, beta0, T,s0,i0,0) )
# print("Epi, Infizierte in D nach", T, " Tagen: ", methods.epidlös(gamma, beta0, T,s0,i0,0) )

# print(EpVerlauf())
# print(PhasenportraitEp())
# print(PhasenportraitEnd())
