#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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