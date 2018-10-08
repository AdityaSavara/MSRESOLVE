# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:07:44 2018

@author: Alex
"""
import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#import the functions from UnitTesterSG
import UnitTesterSG as ut

#imports the KeepOnlySelectedYYYYColumns from XYYYDataFunctionsSG.py
from XYYYDataFunctionsSG import KeepOnlySelectedYYYYColumns

#get the suffix argument for check_results
suffix=ut.returnDigitFromFilename(__file__)

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
def test_Run(allowOverwrite = False):
    #if the user wants to be able to change what the saved outputs are
    if allowOverwrite:
        #This function call is used when this test is run solo as well as by UnitTesterSG
        ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix)
    #this option allows pytest to call the function
    if not allowOverwrite: 
        #this assert statement is required for the pytest module 
        assert ut.check_results(resultObj, resultStr, prefix = '', suffix=suffix, allowOverwrite = False) == True
    
if __name__ == "__main__":
   test_Run(allowOverwrite = True)


