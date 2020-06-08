#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:33:03 2020

@author: be
"""


def ForwardEuler(fun, y0, T,n):
    '''Forwar Euler nimmt als Parameter: \n
        fun, Ableitung der Funktion abhänjgig von x und y \n
        y0,  Anfangswert \n
        T, Endzeitpunkt \n
        n, Anzahl der Iterationsschritte'''
        
    Schritt = y0           
    for i in range(n+1): 
        Schritt = Schritt +  T/n*fun(T/n * i, Schritt)
    return Schritt