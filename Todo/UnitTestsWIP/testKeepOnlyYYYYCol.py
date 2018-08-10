# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:08:08 2018

@author: 3cw
"""
import numpy as np
import XYYYDataFunctionsSG as XYYY

bigArray = np.array([[1, 2, 3, 4, 5],
                     [1, 2, 3, 4, 5],
                     [1, 2, 3, 4, 5],
                     [1, 2, 3, 4, 5]])

headers = ['a','b','c','d','e']
headersToKeep = ['a','e','c']


smallerArray,Headers = XYYY.KeepOnlySelectedYYYYColumns(bigArray,headers,headersToKeep)

print(headers)
print(smallerArray)
