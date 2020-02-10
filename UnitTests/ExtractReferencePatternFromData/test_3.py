# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 13:26:37 2018

@author: Alex
"""
import sys
import os
sys.path.insert(1, os.path.join(os.curdir, os.pardir, "lib"))
sys.path.insert(1, os.path.join(os.curdir, os.pardir))
sys.path.insert(1, os.path.join(os.curdir, os.pardir, os.pardir))
#This test file tests the extractReferencePatternFromData feature

import MSRESOLVE, importlib; importlib.reload(MSRESOLVE)
import UnitTesterSG as ut

#Import the test input file
import test_3_input

#This replaces the globals variables being pointed to in MSRESOLVE to our test_i_input
MSRESOLVE.G = test_3_input

#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)

#The function in main function being tested is extractReferencePatternFromData
MSRESOLVE.main()

#Find the index of the exported data in the MSReference class
index = MSRESOLVE.currentReferenceData.labelToExport.index('ExtractReferencePatternFromData')

#The output will be the dataToExport of the same index as labelToExport
output = MSRESOLVE.currentReferenceData.dataToExport[index]
    
#places the object in a tuple
resultObj = (output)

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
