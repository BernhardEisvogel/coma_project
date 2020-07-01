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

print("             . \n"
   "      .   :   .  Herzlich Willkommen beim CoMa -\n"+
   "  '.   .  :  .   .---       Projekt von Louisa und Bernhard !\n"+
 " ._   '._.-'''-._.'   _._______.-------\n"+
 "   '-..'         '..-' \n"+
" --._ /.==.     .==.\ _.------^          Stay  \n"+
"     ;/_o__\   /_o__\;\n"+
"-----|`#### ) ( ####`|--------.________  Safe and Healthy !!\n"+
"    _:     (\_/)     ;_\n"+
" --'  \   ______     /  '-------\n"+
"   _.-''. ||||||   .''-._  \n"+
"  '    .''-.._..-''.     \  '\n"+
"    .'   '  :  '   '     |. <-- Very angry Virus, because now YOU\n"+
"                                      can predict it's secret moves !\n"+
 "       '    :   '")
def EpVerlauf(epi = True):
    """
    Dieses Programm liest die Datei sirdyn.param und zeichnet für die daraus\n
    ausgelesenen Werte ein Diagramm mit dem Verlauf der Infizierten und An- \n
    fälligen Zahl
    
    Es wird standardmäßig das epidemologische Modell benutzt. Dies kann aber \n
    mit epi = False als Eingabewert umgestellt werden.
    """
    gamma = 0
    beta0 =0
    my =0
    Tgelesen = 0
    s=[]
    i=[]
    s0 = 0
    i0 = 0
    t=[0]

    if(epi == True): 
        gamma, beta0, Tgelesen, s0,i0 = methods.SirLesen()
        s=[s0]
        i=[i0]
        for r in np.arange(1,Tgelesen+0.1,1):
            t.append(r)
            loesung = methods.epidlös(gamma, beta0, 1, s[len(s)-1],i[len(i)-1])
            s.append(loesung[0])
            i.append(loesung[1])

    else:
        gamma, beta0, my, Tgelesen, s0,i0 = methods.SirDynLesen()
        s=[s0]
        i=[i0]
        for r in np.arange(1,Tgelesen+0.1,1):
            t.append(r)
            loesung = methods.endlös(my, gamma, beta0, 1, s[len(s)-1],i[len(i)-1])
                                     
            s.append(loesung[0])
            i.append(loesung[1])
            
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
        arbeitsergebnis = methods.epidlös(gamma,beta0,i,s0,i0)
        e=np.append(e,[str(int((1-arbeitsergebnis[0])*N))])
        r=np.append(r,[str(int((1-arbeitsergebnis[0]-arbeitsergebnis[1])*N))])
    daten=np.vstack((t,e,r))
    daten=daten.T
    np.savetxt("daten.sir.txt",daten,fmt="%20s")

def PhasenportraitEp():
    '''
    Dieses Programm zeichnet ein Phasenportrait des epidemologischen Epidemie-\n
    verlauf. Da

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
            xw,yw= methods.epidlös(gamma, beta0, 3,xw,yw)
        
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
            xw,yw =methods.endlös(my,gamma,beta0,2,xw,yw)
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
    y=[(N-130)/N,130/N]
    for i in range(0,len(inf),1):
        t.append(i)
        I.append(N*(y[1]))
        y=methods.endlös(1/27375,1/6.5,Kontakt[i],i,y[0],y[1])

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
    '''
    Die Funktion erstellt eine Prognose für die Infiziertenzahl am 1.7. \n
    mithilfe des endemischen Modells, ausgehend von der Infiziertenzahl am \n
    27.5. mit der Kontaktrate vom 27.5.
    Returns
    -------
    int, Infiziertenzahl am 1.7.
    '''
    N=83*(10**6)
    my=1/27375
    gamma=1/6.5
    datum=35 # Prognose für den 1.7., 35 Tage nach dem 27.5.
    beta=methods.Kontaktrate(87)[87] #berechnete Kontaktrate vom 27.5.
    s=N-10318-162820 #Anzahl der Anfälligen vom 27.5., berechnet aus den Tabellendaten
    prognose=round(methods.endlös(my,gamma,beta,datum,s/N,10318/N,162820/N)[1]*N)
    print(prognose)


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
    x =[]
    i = []
    s = []
    echteWerte = [1-(math.exp(t))/(math.exp(t) -0.9),(math.exp(t))/(math.exp(t) -0.9),0]
    for n in 1000, 2000, 4000, 8000, 16000:
        approximiert = methods.epidlös(gamma, beta0, t, s0,i0, schritte = n)
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
    x =[]
    i = []
    s = []
    f= []
    echteWerte = [1-(math.exp(t))/(math.exp(t) -0.9),(math.exp(t))/(math.exp(t) -0.9),0]
    for n in 1000, 2000, 4000, 8000, 16000:
        approximiert = methods.epidlös(gamma, beta0, t, s0,i0, schritte = n)
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
    for r in np.arange(0,TBeginnImpfung,1): #Vor der Impfung
        t.append(r)
        loesung = methods.endlös(0, gamma, beta0, 1, s[len(s)-1],i[len(i)-1])
        s.append(loesung[0])
        i.append(loesung[1])
    maximal = 0 # benötigt zur Betrachtung der Belastungsgrenze des Krankenhauses
    for r in np.arange(TBeginnImpfung,TEndeImpfung,1): #während der Impfperiode
        t.append(r)
        loesung, maximalneu = methods.maximalWert(gamma, beta0, 1, s[len(s)-1],
                                                  i[len(i)-1],p)
        if(maximalneu > maximal): maximal = maximalneu
        s.append(loesung[0])
        i.append(loesung[1])
    for r in np.arange(TEndeImpfung,TEndeModell,1): # nach der Impfperiode
        t.append(r)
        loesung, maximalneu = methods.maximalWert(gamma, beta0, 1, s[len(s)-1],
                                                  i[len(i)-1],0)
        if(maximalneu > maximal): maximal = maximalneu
        s.append(loesung[0])
        i.append(loesung[1])
        
    
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
