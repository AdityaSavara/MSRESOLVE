# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:07:44 2018

@author: Alex
"""

#import the functions from UnitTesterSG
from UnitTesterSG import *

#Extracting the digit from the file name to use as prefix/suffix in check_results
def returnDigitFromFilename():
    import os
    filename = os.path.basename(__file__)
    import re
    listOfNumbers = re.findall('\d+',filename)
    extractedDigit = listOfNumbers[0]
    return extractedDigit

#imports the KeepOnlySelectedYYYYColumns from XYYYDataFunctionsSG.py
from XYYYDataFunctionsSG import KeepOnlySelectedYYYYColumns

#get the suffix argument for check_results
suffix = returnDigitFromFilename()

#input
import numpy as np
YYYYData = np.array([[1,2,3,4,5,6,7,8,9,10,11,12,13,14],
                    [2,2,2,2,2,2,2,2,2,22,22,22,22,22],
                    [3,3,3,3,3,3,3,3,3,33,33,33,33,33],
                    [4,4,4,4,4,4,4,4,4,44,44,44,44,44],
                    [5,5,5,5,5,5,5,5,5,55,55,55,55,55]])
headerValues = np.array(['2','18','26','27','28','29','31','39','41','44','45','56','57','70'])
headerValuesToKeep = ['2','18']

#output
output = KeepOnlySelectedYYYYColumns(YYYYData,headerValues,headerValuesToKeep)
resultObj = output

#String is provided
resultStr = str(resultObj)

#run the Unit Tester
check_results(resultObj, resultStr, prefix = '', suffix=suffix)


