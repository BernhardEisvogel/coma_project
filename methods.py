#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:33:03 2020

@author: be
"""


def ForwardEuler(fun, y0, T,n):
    Schritt = y0           
    for i in range(n+1): 
        Schritt = Schritt +  T/n*fun(T/n * i, Schritt)
    return Schritt