# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 16:53:46 2020

@author: Louisa-7440
"""
#Dieses Programm liest die notwendigen Parameter ein

file=open("sir.param","r")
data=[]
for i in file:
    #einzelne Parameter werden getrennt und in Liste gespeichert
    data+=i.split()
    
#konvertieren der Parameter zu floats
for j in (0,len(data)-1):
    data[j]=float(data[j])
    
#Bezeichnung der Parameter
gamma,beta0,my,s0,i0=data[0],data[1],data[2],data[3],data[4]
    
