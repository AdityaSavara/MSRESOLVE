# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 09:23:07 2018

@author: 3cw
"""
import numpy as np

import XYYYDataFunctionsSG

a= np.array([[1, 2, 3, 3, 5],
             [2, 2, 1, 5, 5],
             [3, 2, 1, 7, 5],
             [4, 2, 2, 9, 5],
             [5, 2, 4, 5, 5]])

b= np.array([1,2,3,4])


c= np.array([[1, 2, 1, 3],
             [2, 2, 1, 3],
             [3, 2, 1, 3],
             [4, 2, 1, 3],
             [5, 2, 1, 3]])

d = np.array([1,3,4])



output = XYYYDataFunctionsSG.RemoveSignals(a,b,c,d)

print(output)

